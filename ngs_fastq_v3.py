
from __future__ import annotations
import gzip, json, os, random
from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, List, Literal, Optional
from ngs_prep_v4 import NGSPrepResult, clean_dna, reverse_complement
from utils_bits_v2 import ensure_dir, write_text

@dataclass
class FastqRecord:
    name: str
    sequence: str
    quality: str

@dataclass
class FASTQSimulationConfig:
    coverage:int=20
    paired_end:bool=True
    read_length_nt:int=150
    mean_q:int=34
    sub_rate:float=0.003
    ins_rate:float=0.0002
    del_rate:float=0.0002
    seed:int=7

def _phred_char(q:int)->str:
    return chr(max(2,min(40,int(q)))+33)
def _qual_string(length:int, mean_q:int, rng:random.Random)->str:
    return "".join(_phred_char(int(round(rng.gauss(mean_q,2.0)))) for _ in range(length))
def _mutate_base(base:str, rng:random.Random)->str:
    return rng.choice([b for b in "ACGT" if b != base])

def mutate_read(seq:str, cfg:FASTQSimulationConfig, rng:random.Random):
    seq=clean_dna(seq); out=[]; i=0
    sub=ins=dele=0
    while i < len(seq):
        b=seq[i]
        if rng.random() < cfg.del_rate:
            dele += 1; i += 1; continue
        if rng.random() < cfg.ins_rate:
            out.append(rng.choice("ACGT")); ins += 1
        if rng.random() < cfg.sub_rate:
            out.append(_mutate_base(b, rng)); sub += 1
        else:
            out.append(b)
        i += 1
    return "".join(out), {"sub":sub,"ins":ins,"del":dele, "changed": int((sub+ins+dele)>0)}

def _slice_read(seq:str, read_len:int, mate:Literal["R1","R2"])->str:
    seq=clean_dna(seq)
    if mate=="R1": return seq[:read_len]
    return reverse_complement(seq[-read_len:])

def simulate_fastq_from_library(prep: NGSPrepResult, cfg: Optional[FASTQSimulationConfig]=None)->Dict[str,Any]:
    config=cfg or FASTQSimulationConfig(); rng=random.Random(config.seed)
    r1=[]; r2=[]
    counts={"reads_total":0,"reads_changed":0,"sub":0,"ins":0,"del":0}
    per_fragment=[]
    for frag in prep.fragments:
        frag_changed=0; frag_reads=0; frag_events={"sub":0,"ins":0,"del":0}
        full=frag.full_sequence
        for cov_i in range(config.coverage):
            name_core=f"{frag.fragment_name}|frag={frag.fragment_index}|cov={cov_i+1}"
            seq1, ev1 = mutate_read(_slice_read(full, config.read_length_nt, "R1"), config, rng)
            r1.append(FastqRecord(name=name_core+"/1", sequence=seq1, quality=_qual_string(len(seq1), config.mean_q, rng)))
            frag_reads += 1; frag_changed += ev1["changed"]; 
            for k in ("sub","ins","del"): counts[k] += ev1[k]; frag_events[k]+=ev1[k]
            counts["reads_changed"] += ev1["changed"]; counts["reads_total"] += 1
            if config.paired_end:
                seq2, ev2 = mutate_read(_slice_read(full, config.read_length_nt, "R2"), config, rng)
                r2.append(FastqRecord(name=name_core+"/2", sequence=seq2, quality=_qual_string(len(seq2), config.mean_q, rng)))
                frag_reads += 1; frag_changed += ev2["changed"]; 
                for k in ("sub","ins","del"): counts[k] += ev2[k]; frag_events[k]+=ev2[k]
                counts["reads_changed"] += ev2["changed"]; counts["reads_total"] += 1
        per_fragment.append({
            "fragment_index": int(frag.fragment_index),
            "fragment_name": frag.fragment_name,
            "is_parity_fragment": bool(getattr(frag, "is_parity_fragment", False)),
            "reads_generated": frag_reads,
            "reads_with_errors": frag_changed,
            "error_fraction": (frag_changed/frag_reads if frag_reads else 0.0),
            "substitutions": frag_events["sub"],
            "insertions": frag_events["ins"],
            "deletions": frag_events["del"],
        })
    summary={
        "read_count_r1": len(r1), "read_count_r2": len(r2), "coverage": int(config.coverage), "paired_end": bool(config.paired_end),
        "read_length_nt": int(config.read_length_nt), "mean_q": int(config.mean_q),
        "error_model":{"sub_rate":float(config.sub_rate),"ins_rate":float(config.ins_rate),"del_rate":float(config.del_rate)},
        "error_injection_enabled": bool(config.sub_rate or config.ins_rate or config.del_rate),
        "total_reads": counts["reads_total"], "reads_with_errors": counts["reads_changed"],
        "read_error_fraction": (counts["reads_changed"]/counts["reads_total"] if counts["reads_total"] else 0.0),
        "total_substitutions": counts["sub"], "total_insertions": counts["ins"], "total_deletions": counts["del"],
        "fragment_error_summary": per_fragment,
    }
    return {"r1":r1,"r2":r2,"summary":summary,"config":asdict(config)}

def write_fastq_records(records: Iterable[FastqRecord], path:str, gzip_output:bool=True)->str:
    ensure_dir(os.path.dirname(path))
    text=[]
    for rec in records:
        text += [f"@{rec.name}", rec.sequence, "+", rec.quality]
    payload="\n".join(text)+("\n" if text else "")
    if gzip_output and not path.endswith(".gz"): path += ".gz"
    if path.endswith(".gz"):
        with gzip.open(path,"wt",encoding="utf-8") as f: f.write(payload)
    else:
        write_text(path, payload)
    return path

def export_simulated_fastq(sim_result: Dict[str,Any], out_dir:str, prefix:str="noisy_reads")->Dict[str,str]:
    ensure_dir(out_dir)
    r1=write_fastq_records(sim_result.get("r1",[]), os.path.join(out_dir,f"{prefix}_with_errors_R1.fastq.gz"), True)
    out={"r1_fastq_gz":r1}
    if sim_result.get("r2"):
        r2=write_fastq_records(sim_result.get("r2",[]), os.path.join(out_dir,f"{prefix}_with_errors_R2.fastq.gz"), True)
        out["r2_fastq_gz"]=r2
    meta=os.path.join(out_dir, f"{prefix}_error_summary.json")
    write_text(meta, json.dumps({"summary":sim_result.get("summary",{}),"config":sim_result.get("config",{})}, indent=2, ensure_ascii=False))
    out["summary_json"]=meta
    return out
