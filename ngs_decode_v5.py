
from __future__ import annotations
import gzip, json, os, uuid, zlib
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple
import dna_codec
from compressors_v2 import restore_rep
from ngs_prep_v4 import FragmentRecord, NGSPrepResult, clean_dna, gc_ratio, max_homopolymer, reconstruct_payload_from_fragments, reverse_complement, shannon_entropy, dna_to_bytes, crc_nt, xor_nt_sequences
from utils_bits_v2 import bitstring_to_bytes, ensure_dir, safe_basename, zlib_inflate_until_eof

def load_json(path:str)->Dict[str,Any]:
    with open(path,"r",encoding="utf-8") as f: return json.load(f)
def fragment_from_dict(d:Dict[str,Any])->FragmentRecord:
    return FragmentRecord(**d)
def load_ngs_prep_result_from_json(path:str)->NGSPrepResult:
    p=load_json(path); fr=[fragment_from_dict(x) for x in p.get("fragments",[])]
    return NGSPrepResult(input_sequence_length_nt=int(p.get("input_sequence_length_nt") or 0), config=dict(p.get("config") or {}), summary=dict(p.get("summary") or {}), fragments=fr)

def _open_text_auto(path:str):
    return gzip.open(path,"rt",encoding="utf-8",errors="ignore") if path.endswith(".gz") else open(path,"r",encoding="utf-8",errors="ignore")
def iter_fastq_records(path:str):
    with _open_text_auto(path) as f:
        while True:
            h=f.readline()
            if not h: break
            seq=f.readline(); plus=f.readline(); qual=f.readline()
            if not qual: break
            yield (h.strip()[1:] if h.startswith("@") else h.strip(), clean_dna(seq.strip()), qual.strip())
def parse_header_frag_idx(header:str)->Optional[int]:
    for part in header.split("|"):
        if part.startswith("frag="):
            try: return int(part.split("=",1)[1].split("/",1)[0])
            except: return None
    return None

def dna_stats(seq:str)->Dict[str,Any]:
    seq=clean_dna(seq)
    return {"length_nt":len(seq),"gc_ratio":gc_ratio(seq) if seq else 0.0,"max_homopolymer":max_homopolymer(seq) if seq else 0,"shannon_entropy":shannon_entropy(seq) if seq else 0.0}

def reconstruct_payload_exact(prep:NGSPrepResult)->str:
    return reconstruct_payload_from_fragments(prep.fragments, drop_padding_on_last=True)

def _consensus_base(counter: Counter, fallback:str)->str:
    if not counter: return fallback
    return sorted(counter.items(), key=lambda kv:(-kv[1], kv[0]))[0][0]

def _build_consensus_fragments(prep: NGSPrepResult, fastq_r1_path:str, fastq_r2_path:Optional[str]=None):
    frag_map={int(f.fragment_index):f for f in prep.fragments}
    base_counts: Dict[int, List[Counter]]={}
    total_reads=0
    def ensure_arr(fidx:int, length:int):
        if fidx not in base_counts: base_counts[fidx]=[Counter() for _ in range(length)]
    for header,seq,qual in iter_fastq_records(fastq_r1_path):
        total_reads += 1
        fidx=parse_header_frag_idx(header)
        if fidx is None or fidx not in frag_map: continue
        full=frag_map[fidx].full_sequence; ensure_arr(fidx,len(full))
        n=min(len(seq),len(full))
        for i in range(n): base_counts[fidx][i][seq[i]] += 1
    if fastq_r2_path and os.path.exists(fastq_r2_path):
        for header,seq,qual in iter_fastq_records(fastq_r2_path):
            total_reads += 1
            fidx=parse_header_frag_idx(header)
            if fidx is None or fidx not in frag_map: continue
            full=frag_map[fidx].full_sequence; ensure_arr(fidx,len(full))
            aligned=reverse_complement(seq); n=min(len(aligned),len(full)); start=len(full)-n
            for off in range(n): base_counts[fidx][start+off][aligned[off]] += 1
    rows=[]; frag_payloads={}; frag_status={}
    for fidx in sorted(frag_map):
        frag=frag_map[fidx]
        counts=base_counts.get(fidx, [Counter() for _ in range(len(frag.full_sequence))])
        cons="".join(_consensus_base(counts[i], frag.full_sequence[i]) for i in range(len(frag.full_sequence)))
        left_len=len(frag.forward_primer)+len(frag.sample_barcode)+len(frag.fragment_id)
        right_len=len(frag.checksum)+len(frag.reverse_primer)
        payload=cons[left_len: len(cons)-right_len if right_len else len(cons)]
        expected_len = (frag.payload_end_nt-frag.payload_start_nt) if (frag.payload_end_nt and not frag.is_parity_fragment) else len(frag.payload)
        if expected_len and len(payload)>expected_len: payload=payload[:expected_len]
        support=sum(sum(c.values())>0 for c in counts)
        mismatch_frac = sum(1 for a,b in zip(payload, frag.payload) if a!=b)/max(1,min(len(payload),len(frag.payload))) if frag.payload else 0.0
        crc_expected = crc_nt(payload, len(frag.checksum)) if frag.checksum else ""
        crc_ok = (not frag.checksum) or (crc_expected == frag.checksum)
        frag_payloads[fidx]=payload
        status="Recovered"
        if support==0: status="Missing"
        elif frag.checksum and not crc_ok: status="CRC Failed"
        rows.append({
            "fragment_name": frag.fragment_name,
            "fragment_index": frag.fragment_index,
            "is_parity_fragment": bool(getattr(frag,"is_parity_fragment",False)),
            "supporting_positions": support,
            "difficulty": frag.difficulty,
            "score": frag.score,
            "payload_length_nt": len(frag.payload),
            "recovered": bool(support>0),
            "status": status,
            "mismatch_fraction_vs_manifest": round(float(mismatch_frac),4),
            "crc_ok": crc_ok if frag.checksum else None,
        })
        frag_status[fidx]=status
    return frag_payloads, rows, total_reads

def recover_payload_from_fastq(prep:NGSPrepResult, fastq_r1_path:str, fastq_r2_path:Optional[str]=None, method:str="Consensus reconstruction")->Tuple[str,Dict[str,Any]]:
    payloads, rows, total_reads = _build_consensus_fragments(prep, fastq_r1_path, fastq_r2_path)
    protection_mode = str((prep.config or {}).get("protection_mode") or "None")
    rescued=0
    if protection_mode == "Parity Fragments" and method in ("Consensus + checksum validation","ECC-assisted recovery","Consensus + manifest-guided rescue"):
        # Single-erasure correction per parity group: if exactly one data fragment missing/CRC-failed.
        data_frags=[f for f in prep.fragments if not getattr(f,"is_parity_fragment",False)]
        parity_frags=[f for f in prep.fragments if getattr(f,"is_parity_fragment",False)]
        for parity in parity_frags:
            covers=[int(x) for x in (parity.parity_covers or [])]
            bad=[idx for idx in covers if idx not in payloads or rows[[r["fragment_index"] for r in rows].index(idx)]["status"] in ("Missing","CRC Failed")]
            if len(bad)==1:
                miss=bad[0]
                others=[payloads[idx] for idx in covers if idx != miss and idx in payloads]
                parity_payload=payloads.get(parity.fragment_index, parity.payload)
                if others and parity_payload:
                    rec = xor_nt_sequences(others+[parity_payload])
                    expected_len=len(next(f.payload for f in prep.fragments if f.fragment_index==miss))
                    rec=rec[:expected_len]
                    payloads[miss]=rec
                    for row in rows:
                        if row["fragment_index"]==miss:
                            row["status"]="Rescued by parity"
                            row["recovered"]=True
                            rescued += 1
    # Build final payload from non-parity fragments
    ordered=[f for f in sorted(prep.fragments, key=lambda x:x.fragment_index) if not getattr(f,"is_parity_fragment",False)]
    parts=[]
    for frag in ordered:
        payload=clean_dna(payloads.get(frag.fragment_index) or "")
        expected_len=max(0,int(frag.payload_end_nt)-int(frag.payload_start_nt)) if frag.payload_end_nt else len(frag.payload)
        if expected_len and len(payload)>expected_len: payload=payload[:expected_len]
        parts.append(payload)
    payload="".join(parts)
    recovered=sum(1 for r in rows if r["recovered"] and not r.get("is_parity_fragment"))
    data_total=sum(1 for r in rows if not r.get("is_parity_fragment"))
    report={
        "mode":"recovery",
        "protection_mode": protection_mode,
        "total_reads": int(total_reads),
        "total_fragments": data_total,
        "recovered_fragments": recovered,
        "recovery_fraction": (recovered/data_total if data_total else 0.0),
        "rescued_fragments": int(rescued),
        "rows": rows,
    }
    return payload, report

def precheck_decode_payload(dna_payload:str, *, scheme_name:str, mode_codec:str, seed:str, init_dimer:str, whiten:bool, remove_leading_one:bool, target_gc:float, w_gc:float, w_motif:float, ks:Tuple[int,...])->Dict[str,Any]:
    dna=clean_dna(dna_payload)
    if not dna: return {"ok":False,"error":"Empty DNA payload"}
    try:
        bits,digits = dna_codec.decode_dna_to_bits(dna, scheme_name=scheme_name, mode=mode_codec, seed=seed, init_dimer=init_dimer, remove_leading_one=bool(remove_leading_one), whiten=bool(whiten), target_gc=float(target_gc), w_gc=float(w_gc), w_motif=float(w_motif), ks=ks)
        decoded_buf,pad_bits=bitstring_to_bytes(bits, pad_to_byte=True)
        inner,infl=zlib_inflate_until_eof(decoded_buf)
        z_ok=bool(infl.get("eof")) and (infl.get("error") is None)
        return {"ok":bool(z_ok),"decoded_bits_len":len(bits),"decoded_bytes_len":len(decoded_buf),"inner_bytes_len":len(inner),"pad_bits_to_byte":pad_bits,"digits_len":(len(digits) if isinstance(digits,list) else None),"zlib":{"eof":bool(infl.get("eof")),"unused_tail_len_bytes":int(infl.get("unused_tail_len",0) or 0),"error":infl.get("error"),"integrity_ok":bool(z_ok)}}
    except Exception as e:
        return {"ok":False,"error":str(e)}

def decode_dna_payload_to_file(dna_payload:str, *, scheme_name:str, mode_codec:str, seed:str, init_dimer:str, whiten:bool, remove_leading_one:bool, target_gc:float, w_gc:float, w_motif:float, ks:Tuple[int,...], preferred_stem:str="recovered_from_fastq")->Tuple[Dict[str,Any], Optional[str]]:
    pre=precheck_decode_payload(dna_payload, scheme_name=scheme_name, mode_codec=mode_codec, seed=seed, init_dimer=init_dimer, whiten=whiten, remove_leading_one=remove_leading_one, target_gc=target_gc, w_gc=w_gc, w_motif=w_motif, ks=ks)
    if not pre.get("ok"): return pre, None
    dna=clean_dna(dna_payload)
    bits,digits = dna_codec.decode_dna_to_bits(dna, scheme_name=scheme_name, mode=mode_codec, seed=seed, init_dimer=init_dimer, remove_leading_one=bool(remove_leading_one), whiten=bool(whiten), target_gc=float(target_gc), w_gc=float(w_gc), w_motif=float(w_motif), ks=ks)
    decoded_buf,pad_bits=bitstring_to_bytes(bits, pad_to_byte=True)
    inner,infl=zlib_inflate_until_eof(decoded_buf)
    z_ok=bool(infl.get("eof")) and (infl.get("error") is None)
    out_dir=os.path.join("recovery_out","ngs_decode",str(uuid.uuid4())); ensure_dir(out_dir)
    restored_file=None; restore_meta=None
    if z_ok:
        restored_file, restore_meta = restore_rep(inner, out_dir=out_dir, preferred_stem=safe_basename(preferred_stem, fallback="recovered_from_fastq"))
    stats={"ok":bool(z_ok),"dna_nt":len(dna),"decoded_bits_len":len(bits),"decoded_bytes_len":len(decoded_buf),"pad_bits_to_byte":pad_bits,"digits_len":(len(digits) if isinstance(digits,list) else None),"zlib":{"eof":bool(infl.get('eof')),"unused_tail_len_bytes":int(infl.get('unused_tail_len',0) or 0),"error":infl.get('error'),"integrity_ok":z_ok},"inner_bytes_len":len(inner),"restored_file":restored_file,"restore_meta":restore_meta}
    return stats, restored_file
