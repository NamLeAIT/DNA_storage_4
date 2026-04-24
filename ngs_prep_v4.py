
from __future__ import annotations

import csv, hashlib, json, math, os, zlib
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple

from utils_bits_v2 import ensure_dir, write_text

DNAAlphabet = set("ACGT")
MAP = "ACGT"
INV = {b:i for i,b in enumerate(MAP)}

def clean_dna(seq: str) -> str:
    if seq is None:
        return ""
    return "".join(ch for ch in str(seq).upper() if ch in DNAAlphabet)

def reverse_complement(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGT","TGCA"))[::-1]

def gc_ratio(seq: str) -> float:
    seq = clean_dna(seq); return (sum(ch in "GC" for ch in seq)/len(seq)) if seq else 0.0

def shannon_entropy(seq: str) -> float:
    seq = clean_dna(seq)
    if not seq: return 0.0
    n=len(seq); ent=0.0
    for b in "ACGT":
        c=seq.count(b)
        if c:
            p=c/n; ent -= p*math.log2(p)
    return ent

def max_homopolymer(seq: str) -> int:
    seq = clean_dna(seq)
    if not seq: return 0
    cur=mx=1
    for i in range(1,len(seq)):
        if seq[i]==seq[i-1]:
            cur += 1; mx=max(mx,cur)
        else:
            cur=1
    return mx

def homopolymer_count(seq: str, min_len: int=2) -> int:
    seq=clean_dna(seq)
    if not seq: return 0
    cur=1; cnt=0
    for i in range(1,len(seq)):
        if seq[i]==seq[i-1]:
            cur += 1
        else:
            if cur >= min_len: cnt += 1
            cur = 1
    if cur >= min_len: cnt += 1
    return cnt

def has_forbidden_motif(seq: str, motifs: Tuple[str,...]) -> bool:
    seq=clean_dna(seq)
    return any(m in seq for m in motifs)

def simple_self_complementarity_score(seq: str, k: int=6) -> int:
    seq=clean_dna(seq)
    if len(seq) < k: return 0
    kmers=set(seq[i:i+k] for i in range(len(seq)-k+1))
    return sum(1 for x in kmers if reverse_complement(x) in kmers)

def _component_score(seq: str, *, target_gc: float=0.50, hp_limit: int=3, forbidden_motifs: Tuple[str,...]=()) -> float:
    seq=clean_dna(seq)
    if not seq: return 100.0
    score=100.0
    score -= abs(gc_ratio(seq)-target_gc)*60.0
    hp=max_homopolymer(seq)
    if hp > hp_limit:
        score -= min(35.0,(hp-hp_limit)*8.0)
    ent=shannon_entropy(seq)
    if ent < 1.8:
        score -= min(20.0,(1.8-ent)*40.0)
    score -= min(12.0,simple_self_complementarity_score(seq, k=min(6,max(3,len(seq)//3))))
    if has_forbidden_motif(seq, forbidden_motifs):
        score -= 12.0
    return max(0.0,min(100.0,score))

def constrained_nt_string(length_nt: int, *, seed_int: int, left_context: str="", right_context: str="", target_gc: float=0.50, forbidden_motifs: Tuple[str,...]=("AAAAAA","CCCCCC","GGGGGG","TTTTTT"), soft_hp_limit: int=2) -> str:
    if length_nt <= 0: return ""
    h = hashlib.sha256(f"{seed_int}|{left_context}|{right_context}|{length_nt}".encode()).digest()
    n = int.from_bytes(h, "big")
    preferred = [MAP[(n >> (2*i)) & 0b11] for i in range(length_nt)]
    out=[]; current=clean_dna(left_context)
    for i in range(length_nt):
        candidates=[preferred[i]]+[b for b in MAP if b != preferred[i]]
        best_b, best_s=None, None
        for b in candidates:
            trial=current+"".join(out)+b
            if i == length_nt-1 and right_context:
                trial += clean_dna(right_context)[:1]
            s = 100.0 - _component_score(trial, target_gc=target_gc, hp_limit=soft_hp_limit, forbidden_motifs=forbidden_motifs)
            if best_s is None or s < best_s:
                best_b, best_s = b, s
        out.append(best_b or "A")
    return "".join(out)

def safe_padding(length_nt: int, *, seed_seq: str="", left_context: str="", right_context: str="", target_gc: float=0.50, forbidden_motifs: Tuple[str,...]=("AAAAAA","CCCCCC","GGGGGG","TTTTTT")) -> str:
    return constrained_nt_string(length_nt, seed_int=int(hashlib.sha256((seed_seq+left_context+right_context).encode()).hexdigest()[:12],16), left_context=left_context, right_context=right_context, target_gc=target_gc, forbidden_motifs=forbidden_motifs, soft_hp_limit=2)

def dna_to_bytes(seq: str) -> bytes:
    seq=clean_dna(seq)
    bits=[]
    for b in seq:
        v=INV[b]
        bits.append((v>>1)&1); bits.append(v&1)
    out=bytearray()
    for i in range(0,len(bits),8):
        chunk=bits[i:i+8]
        while len(chunk)<8: chunk.append(0)
        val=0
        for bit in chunk: val=(val<<1)|bit
        out.append(val)
    return bytes(out)

def bytes_to_dna(data: bytes, nt_len: int) -> str:
    bits=[]
    for b in data:
        for i in range(7,-1,-1):
            bits.append((b>>i)&1)
    out=[]
    for i in range(0, min(len(bits), nt_len*2), 2):
        v=(bits[i]<<1) | (bits[i+1] if i+1 < len(bits) else 0)
        out.append(MAP[v])
    while len(out) < nt_len:
        out.append("A")
    return "".join(out[:nt_len])

def crc_nt(payload_nt: str, crc_nt_len: int) -> str:
    if crc_nt_len <= 0: return ""
    crc = zlib.crc32(dna_to_bytes(payload_nt))
    byte_len=max(1, math.ceil(crc_nt_len/4))
    data=crc.to_bytes(4,'big')[-byte_len:]
    return bytes_to_dna(data, crc_nt_len)

def xor_nt_sequences(seqs: List[str]) -> str:
    if not seqs: return ""
    L=max(len(s) for s in seqs)
    out=[]
    for i in range(L):
        x=0
        for s in seqs:
            b = s[i] if i < len(s) else "A"
            x ^= INV.get(b,0)
        out.append(MAP[x])
    return "".join(out)

@dataclass
class PrimerConfig:
    forward: str = "ACACGACGCTCTTCCGATCT"
    reverse: str = "AGATCGGAAGAGCACACGTCT"
    def cleaned(self)->"PrimerConfig":
        return PrimerConfig(clean_dna(self.forward), clean_dna(self.reverse))

@dataclass
class NGSPrepConfig:
    payload_length_nt: int = 80
    fragment_id_length_nt: int = 8
    sample_barcode: str = ""
    checksum_length_nt: int = 0
    primer: PrimerConfig = field(default_factory=PrimerConfig)
    last_fragment_policy: Literal["keep_shorter","pad_to_full"] = "keep_shorter"
    target_gc: float = 0.50
    gc_soft_range: Tuple[float,float] = (0.40,0.60)
    homopolymer_soft_limit: int = 4
    forbidden_motifs: Tuple[str,...] = ("AAAAAA","CCCCCC","GGGGGG","TTTTTT")
    protection_mode: Literal["None","CRC","Parity Fragments"] = "None"
    parity_group_size: int = 4

    def normalized(self)->"NGSPrepConfig":
        p=self.primer.cleaned()
        return NGSPrepConfig(
            payload_length_nt=max(1,int(self.payload_length_nt)),
            fragment_id_length_nt=max(0,int(self.fragment_id_length_nt)),
            sample_barcode=clean_dna(self.sample_barcode),
            checksum_length_nt=max(0,int(self.checksum_length_nt)),
            primer=p,
            last_fragment_policy=self.last_fragment_policy,
            target_gc=float(self.target_gc),
            gc_soft_range=(float(self.gc_soft_range[0]), float(self.gc_soft_range[1])),
            homopolymer_soft_limit=int(self.homopolymer_soft_limit),
            forbidden_motifs=tuple(clean_dna(x) for x in self.forbidden_motifs if clean_dna(x)),
            protection_mode=str(self.protection_mode),
            parity_group_size=max(2,int(self.parity_group_size)),
        )

@dataclass
class FragmentRecord:
    fragment_name: str
    fragment_index: int
    total_fragments: int
    full_sequence: str
    forward_primer: str
    sample_barcode: str
    fragment_id: str
    payload: str
    checksum: str
    reverse_primer: str
    payload_length_nt: int
    fragment_length_nt: int
    is_padded: bool
    is_last_fragment: bool
    gc_ratio: float
    max_homopolymer: int
    homopolymer_count: int
    shannon_entropy: float
    self_complementarity_score: int
    forbidden_motif_hit: bool
    score: int
    difficulty: str
    payload_start_nt: int
    payload_end_nt: int
    payload_gc_ratio: float = 0.0
    payload_max_homopolymer: int = 0
    payload_shannon_entropy: float = 0.0
    protection_mode: str = "None"
    is_parity_fragment: bool = False
    parity_group_index: int = 0
    parity_covers: List[int] = field(default_factory=list)
    crc_valid_expected: bool = False

@dataclass
class NGSPrepResult:
    input_sequence_length_nt: int
    config: Dict[str,Any]
    summary: Dict[str,Any]
    fragments: List[FragmentRecord]

def _score_fragment(seq: str, *, target_gc: float, gc_soft_range: Tuple[float,float], hp_limit: int, forbidden_motifs: Tuple[str,...]) -> Dict[str,Any]:
    gc=gc_ratio(seq); hp=max_homopolymer(seq); hp_cnt=homopolymer_count(seq,2); ent=shannon_entropy(seq); selfc=simple_self_complementarity_score(seq,6); motif=has_forbidden_motif(seq,forbidden_motifs)
    score=100.0
    lo,hi=gc_soft_range
    if gc<lo: score -= min(20.0,((lo-gc)/max(lo,1e-9))*20.0)
    elif gc>hi: score -= min(20.0,((gc-hi)/max(1-hi,1e-9))*20.0)
    else: score -= abs(gc-target_gc)*10.0
    if hp > hp_limit: score -= min(25.0,(hp-hp_limit)*6.0)
    score -= min(10.0,hp_cnt*0.25)
    if ent<1.8: score -= min(20.0,(1.8-ent)*40.0)
    score -= min(10.0,selfc*1.0)
    if motif: score -= 10.0
    score_i=max(0,min(100,int(round(score))))
    difficulty="Easy" if score_i>=85 else ("Medium" if score_i>=65 else "Hard")
    return dict(gc_ratio=gc,max_homopolymer=hp,homopolymer_count=hp_cnt,shannon_entropy=ent,self_complementarity_score=selfc,forbidden_motif_hit=motif,score=score_i,difficulty=difficulty)

def _make_fragment(cfg: NGSPrepConfig, *, idx: int, total: int, payload: str, payload_start: int, payload_end: int, is_last: bool, is_padded: bool, is_parity: bool=False, parity_group_index: int=0, parity_covers: Optional[List[int]]=None) -> FragmentRecord:
    frag_id = constrained_nt_string(cfg.fragment_id_length_nt, seed_int=(100000+idx if is_parity else idx), left_context=cfg.primer.forward+cfg.sample_barcode, right_context=payload[:4], target_gc=cfg.target_gc, forbidden_motifs=cfg.forbidden_motifs, soft_hp_limit=2) if cfg.fragment_id_length_nt>0 else ""
    checksum = ""
    if cfg.protection_mode == "CRC":
        crc_len = max(4, cfg.checksum_length_nt or 8)
        checksum = crc_nt(payload, crc_len)
    full = cfg.primer.forward + cfg.sample_barcode + frag_id + payload + checksum + cfg.primer.reverse
    m=_score_fragment(full, target_gc=cfg.target_gc, gc_soft_range=cfg.gc_soft_range, hp_limit=cfg.homopolymer_soft_limit, forbidden_motifs=cfg.forbidden_motifs)
    return FragmentRecord(
        fragment_name=f"seq_{idx}",
        fragment_index=idx,
        total_fragments=total,
        full_sequence=full,
        forward_primer=cfg.primer.forward,
        sample_barcode=cfg.sample_barcode,
        fragment_id=frag_id,
        payload=payload,
        checksum=checksum,
        reverse_primer=cfg.primer.reverse,
        payload_length_nt=len(payload),
        fragment_length_nt=len(full),
        is_padded=is_padded,
        is_last_fragment=is_last,
        payload_start_nt=payload_start,
        payload_end_nt=payload_end,
        payload_gc_ratio=gc_ratio(payload),
        payload_max_homopolymer=max_homopolymer(payload),
        payload_shannon_entropy=shannon_entropy(payload),
        protection_mode=cfg.protection_mode,
        is_parity_fragment=is_parity,
        parity_group_index=parity_group_index,
        parity_covers=list(parity_covers or []),
        crc_valid_expected=(cfg.protection_mode=="CRC" and not is_parity),
        **m
    )

def prepare_ngs_library(dna_sequence: str, config: Optional[NGSPrepConfig]=None) -> NGSPrepResult:
    seq=clean_dna(dna_sequence)
    cfg=(config or NGSPrepConfig()).normalized()
    if not seq: raise ValueError("dna_sequence is empty after cleaning")
    payload_len=int(cfg.payload_length_nt)
    parts=[seq[i:i+payload_len] for i in range(0,len(seq),payload_len)]
    rows=[]
    for idx,payload0 in enumerate(parts, start=1):
        payload=payload0; is_last=idx==len(parts); is_padded=False
        if len(payload) < payload_len and cfg.last_fragment_policy == "pad_to_full":
            payload = payload + safe_padding(payload_len-len(payload), seed_seq=payload, left_context=payload, right_context=cfg.primer.reverse, target_gc=cfg.target_gc, forbidden_motifs=cfg.forbidden_motifs)
            is_padded=True
        rows.append(_make_fragment(cfg, idx=idx, total=0, payload=payload, payload_start=(idx-1)*payload_len, payload_end=min(idx*payload_len,len(seq)), is_last=is_last, is_padded=is_padded))
    # block parity over data payloads
    if cfg.protection_mode == "Parity Fragments":
        group_size=max(2,int(cfg.parity_group_size))
        next_idx=len(rows)+1
        for gstart in range(0, len(rows), group_size):
            group=rows[gstart:gstart+group_size]
            if not group: continue
            parity_payload = xor_nt_sequences([r.payload for r in group])
            rec=_make_fragment(cfg, idx=next_idx, total=0, payload=parity_payload, payload_start=0, payload_end=0, is_last=False, is_padded=False, is_parity=True, parity_group_index=(gstart//group_size)+1, parity_covers=[r.fragment_index for r in group])
            rows.append(rec); next_idx += 1
    total=len(rows)
    for r in rows: r.total_fragments=total
    easy=sum(1 for r in rows if r.difficulty=="Easy"); medium=sum(1 for r in rows if r.difficulty=="Medium"); hard=sum(1 for r in rows if r.difficulty=="Hard")
    total_library_nt=sum(r.fragment_length_nt for r in rows)
    data_count=sum(1 for r in rows if not r.is_parity_fragment)
    parity_count=sum(1 for r in rows if r.is_parity_fragment)
    summary=dict(
        total_fragments=total,
        data_fragments=data_count,
        parity_fragments=parity_count,
        full_fragments=len(seq)//payload_len,
        remainder_nt=len(seq)%payload_len,
        padded_last_fragment=bool(rows and rows[data_count-1].is_padded),
        standard_fragment_length_nt=(rows[0].fragment_length_nt if rows else 0),
        last_fragment_length_nt=(rows[data_count-1].fragment_length_nt if rows else 0),
        avg_score=(sum(r.score for r in rows)/len(rows) if rows else 0.0),
        difficulty_distribution={"easy":easy,"medium":medium,"hard":hard},
        total_library_nt=total_library_nt,
        expansion_factor=(total_library_nt/len(seq) if len(seq) else None),
    )
    return NGSPrepResult(input_sequence_length_nt=len(seq), config=asdict(cfg), summary=summary, fragments=rows)

def export_ngs_library(result: NGSPrepResult, out_dir: str, prefix: str="ngs_library") -> Dict[str,str]:
    ensure_dir(out_dir)
    fasta_path=os.path.join(out_dir,f"{prefix}.fasta"); csv_path=os.path.join(out_dir,f"{prefix}.csv"); json_path=os.path.join(out_dir,f"{prefix}.json")
    fasta=[]
    for rec in result.fragments:
        title = rec.fragment_name if not rec.is_parity_fragment else f"{rec.fragment_name}|parity_group={rec.parity_group_index}"
        fasta += [f">{title}", rec.full_sequence]
    write_text(fasta_path, "\n".join(fasta)+"\n")
    fieldnames=list(asdict(result.fragments[0]).keys()) if result.fragments else []
    with open(csv_path,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
        for r in result.fragments: w.writerow(asdict(r))
    payload={"input_sequence_length_nt":result.input_sequence_length_nt,"config":result.config,"summary":result.summary,"fragments":[asdict(x) for x in result.fragments]}
    write_text(json_path, json.dumps(payload, indent=2, ensure_ascii=False))
    return {"fasta":fasta_path,"csv":csv_path,"json":json_path}

def reconstruct_payload_from_fragments(fragments: List[FragmentRecord], *, drop_padding_on_last: bool=True) -> str:
    ordered=sorted([x for x in fragments if not x.is_parity_fragment], key=lambda x:x.fragment_index)
    parts=[]
    for rec in ordered:
        payload=clean_dna(rec.payload)
        expected=max(0,int(rec.payload_end_nt)-int(rec.payload_start_nt))
        if expected and len(payload)>expected: payload=payload[:expected]
        parts.append(payload)
    return "".join(parts)
