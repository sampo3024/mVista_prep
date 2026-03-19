#!/usr/bin/env python3
"""mVISTA input helper GUI.

Builds ortholog-centered regions (promoter or gene+flanks) from genome FASTA + GFF,
and exports files that can be uploaded to mVISTA web.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Dict, List, Optional, Set, Tuple


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class SpeciesInput:
    species: str
    fasta_path: Path
    gff_path: Path
    gene_id: str


@dataclass
class GeneFeature:
    seqid: str
    start: int
    end: int
    strand: str
    raw_attrs: Dict[str, str]


@dataclass
class RegionContext:
    seqid: str
    region_start: int
    region_end: int
    extract_strand: str
    output_seqid: str


@dataclass
class GffRecord:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: Dict[str, str]


class FastaDB:
    def __init__(self, path: Path):
        self.path = path
        self.records: Dict[str, str] = {}
        self._load()

    def _load(self) -> None:
        current = None
        chunks: List[str] = []
        with self.path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current is not None:
                        self.records[current] = "".join(chunks)
                    current = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if current is not None:
                self.records[current] = "".join(chunks)

    def get_subseq(self, seqid: str, start: int, end: int, strand: str) -> str:
        if seqid not in self.records:
            raise KeyError(f"Sequence ID '{seqid}' not found in {self.path}")
        seq = self.records[seqid]
        left = max(1, start)
        right = min(len(seq), end)
        if left > right:
            return ""
        out = seq[left - 1 : right]
        if strand == "-":
            out = out.translate(DNA_COMP)[::-1]
        return out


ATTR_DELIMS = re.compile(r"\s*;\s*")


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for pair in ATTR_DELIMS.split(attr_text.strip(";")):
        if not pair:
            continue
        if "=" in pair:
            key, value = pair.split("=", 1)
        elif " " in pair:
            key, value = pair.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def match_gene(attrs: Dict[str, str], target_gene: str) -> bool:
    target = target_gene.strip()
    keys = ["ID", "Name", "gene", "gene_id", "locus_tag", "transcript_id"]
    for key in keys:
        value = attrs.get(key)
        if value and value == target:
            return True
    return False


def find_gene_feature(gff_path: Path, gene_id: str) -> GeneFeature:
    best: Optional[GeneFeature] = None
    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_txt = fields
            attrs = parse_attrs(attrs_txt)
            if not match_gene(attrs, gene_id):
                continue
            if feature_type not in {"gene", "mRNA", "transcript", "CDS", "exon"}:
                continue
            feat = GeneFeature(
                seqid=seqid,
                start=int(start),
                end=int(end),
                strand=strand,
                raw_attrs=attrs,
            )
            # Prefer "gene" over other feature types.
            if feature_type == "gene":
                return feat
            if best is None:
                best = feat
    if best is not None:
        return best
    raise ValueError(f"Gene '{gene_id}' not found in {gff_path}")


def region_for_feature(feature: GeneFeature, mode: str, upstream: int, downstream: int) -> Tuple[int, int]:
    if mode == "Promoter (upstream from TSS)":
        if feature.strand == "-":
            start = feature.end + 1
            end = feature.end + upstream
        else:
            start = feature.start - upstream
            end = feature.start - 1
    else:
        # Gene body + flanks.
        if feature.strand == "-":
            start = feature.start - downstream
            end = feature.end + upstream
        else:
            start = feature.start - upstream
            end = feature.end + downstream
    return (max(1, start), max(1, end))


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def safe_seqid(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> Optional[Tuple[int, int]]:
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    if s > e:
        return None
    return (s, e)


def remap_to_region(start: int, end: int, ctx: RegionContext) -> Tuple[int, int]:
    if ctx.extract_strand == "-":
        mapped_start = ctx.region_end - end + 1
        mapped_end = ctx.region_end - start + 1
    else:
        mapped_start = start - ctx.region_start + 1
        mapped_end = end - ctx.region_start + 1
    return (mapped_start, mapped_end)


def remap_strand(original: str, extract_strand: str) -> str:
    if original not in {"+", "-"}:
        return original
    if extract_strand == "-":
        return "+" if original == "-" else "-"
    return original


def parse_gff_records(gff_path: Path, seqid_filter: str) -> List[GffRecord]:
    records: List[GffRecord] = []
    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs_txt = fields
            if seqid != seqid_filter:
                continue
            try:
                feat_start = int(start)
                feat_end = int(end)
            except ValueError:
                continue
            records.append(
                GffRecord(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start=feat_start,
                    end=feat_end,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attrs=parse_attrs(attrs_txt),
                )
            )
    return records


def normalize_child_type(feature_type: str) -> Optional[str]:
    ft = feature_type.lower()
    if ft == "exon":
        return "exon"
    if ft in {"utr", "five_prime_utr", "three_prime_utr"}:
        return "utr"
    return None


def record_name(rec: GffRecord) -> str:
    return rec.attrs.get("Name") or rec.attrs.get("gene") or rec.attrs.get("ID") or "gene"


def infer_group_id(rec: GffRecord) -> Optional[str]:
    parents = rec.attrs.get("Parent", "")
    parent_first = ""
    if parents:
        parent_first = [x.strip() for x in parents.split(",") if x.strip()][0]
    return (
        rec.attrs.get("gene_id")
        or rec.attrs.get("gene")
        or rec.attrs.get("locus_tag")
        or parent_first
        or rec.attrs.get("transcript_id")
        or rec.attrs.get("ID")
    )


def gene_id_candidates(rec: GffRecord) -> Set[str]:
    keys = ["ID", "Name", "gene", "gene_id", "locus_tag", "transcript_id"]
    return {rec.attrs[k] for k in keys if rec.attrs.get(k)}


def child_belongs_to_gene(child: GffRecord, gene_ids: Set[str]) -> bool:
    parents = {x.strip() for x in child.attrs.get("Parent", "").split(",") if x.strip()}
    if parents & gene_ids:
        return True
    for key in ["gene_id", "gene", "locus_tag", "transcript_id", "ID", "Name"]:
        v = child.attrs.get(key)
        if v and v in gene_ids:
            return True
    return False


def count_region_records(records: List[GffRecord], ctx: RegionContext) -> Tuple[int, int]:
    total = 0
    child_like = 0
    for rec in records:
        if overlap(rec.start, rec.end, ctx.region_start, ctx.region_end) is None:
            continue
        total += 1
        if normalize_child_type(rec.feature_type) is not None:
            child_like += 1
    return total, child_like


def export_remapped_annotations(
    gff_path: Path,
    ctx: RegionContext,
    feature_types: Set[str],
    gff_out: Path,
    bed_out: Path,
) -> int:
    out_gff_lines: List[str] = ["##gff-version 3"]
    out_bed_lines: List[str] = []
    kept = 0

    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs_txt = fields
            if seqid != ctx.seqid:
                continue
            if feature_types and feature_type.lower() not in feature_types:
                continue

            try:
                feat_start = int(start)
                feat_end = int(end)
            except ValueError:
                continue

            ov = overlap(feat_start, feat_end, ctx.region_start, ctx.region_end)
            if ov is None:
                continue

            clip_start, clip_end = ov
            mapped_start, mapped_end = remap_to_region(clip_start, clip_end, ctx)
            mapped_strand = remap_strand(strand, ctx.extract_strand)
            attrs = parse_attrs(attrs_txt)
            attr_pairs = [f"{k}={v}" for k, v in attrs.items()]
            if "orig_coords" not in attrs:
                attr_pairs.append(f"orig_coords={ctx.seqid}:{clip_start}-{clip_end}")
            attrs_out = ";".join(attr_pairs) if attr_pairs else "."

            out_gff_lines.append(
                "\t".join(
                    [
                        ctx.output_seqid,
                        source,
                        feature_type,
                        str(mapped_start),
                        str(mapped_end),
                        score,
                        mapped_strand,
                        phase,
                        attrs_out,
                    ]
                )
            )

            bed_start = mapped_start - 1
            bed_end = mapped_end
            name = attrs.get("Name") or attrs.get("ID") or feature_type
            out_bed_lines.append(
                "\t".join(
                    [
                        ctx.output_seqid,
                        str(bed_start),
                        str(bed_end),
                        name,
                        "0",
                        mapped_strand if mapped_strand in {"+", "-"} else ".",
                    ]
                )
            )
            kept += 1

    gff_out.write_text("\n".join(out_gff_lines) + "\n")
    bed_out.write_text("\n".join(out_bed_lines) + ("\n" if out_bed_lines else ""))
    return kept


def export_mvista_simple_annotation(gff_path: Path, ctx: RegionContext, out_path: Path) -> int:
    records = parse_gff_records(gff_path, ctx.seqid)
    gene_like = {"gene", "mrna", "transcript"}
    in_region: List[GffRecord] = []

    for rec in records:
        if overlap(rec.start, rec.end, ctx.region_start, ctx.region_end) is not None:
            in_region.append(rec)

    lines: List[str] = []
    genes_written = 0
    gene_records = [r for r in in_region if r.feature_type.lower() in gene_like]
    gene_records.sort(key=lambda r: (r.start, r.end))
    child_records = [r for r in in_region if normalize_child_type(r.feature_type) is not None]

    if gene_records:
        for gene in gene_records:
            ov = overlap(gene.start, gene.end, ctx.region_start, ctx.region_end)
            if ov is None:
                continue
            clip_start, clip_end = ov
            gstart, gend = remap_to_region(clip_start, clip_end, ctx)
            gstrand = remap_strand(gene.strand, ctx.extract_strand)
            sign = ">" if gstrand != "-" else "<"
            lines.append(f"{sign} {gstart} {gend} {record_name(gene)}")
            genes_written += 1

            ids = gene_id_candidates(gene)
            child_entries: List[Tuple[int, int, str]] = []
            for child in child_records:
                if not child_belongs_to_gene(child, ids):
                    continue
                cov = overlap(child.start, child.end, ctx.region_start, ctx.region_end)
                if cov is None:
                    continue
                ctype = normalize_child_type(child.feature_type)
                if ctype is None:
                    continue
                cstart, cend = remap_to_region(cov[0], cov[1], ctx)
                child_entries.append((cstart, cend, ctype))

            child_entries.sort(key=lambda x: (x[0], x[1], x[2]))
            for cstart, cend, ctype in child_entries:
                lines.append(f"{cstart} {cend} {ctype}")
            lines.append("")
    else:
        # Fallback for GTF-like data with no explicit gene/transcript rows:
        # group exon/UTR records by gene_id/locus_tag/Parent and emit pseudo-gene blocks.
        grouped: Dict[str, List[GffRecord]] = {}
        for child in child_records:
            gid = infer_group_id(child)
            if gid:
                grouped.setdefault(gid, []).append(child)

        ordered_groups = sorted(grouped.items(), key=lambda kv: min(r.start for r in kv[1]))
        for gid, group in ordered_groups:
            min_s = min(r.start for r in group)
            max_e = max(r.end for r in group)
            gov = overlap(min_s, max_e, ctx.region_start, ctx.region_end)
            if gov is None:
                continue
            gstart, gend = remap_to_region(gov[0], gov[1], ctx)
            strand0 = group[0].strand if group else "."
            gstrand = remap_strand(strand0, ctx.extract_strand)
            sign = ">" if gstrand != "-" else "<"
            lines.append(f"{sign} {gstart} {gend} {gid}")
            genes_written += 1

            child_entries: List[Tuple[int, int, str]] = []
            for child in group:
                cov = overlap(child.start, child.end, ctx.region_start, ctx.region_end)
                if cov is None:
                    continue
                ctype = normalize_child_type(child.feature_type)
                if ctype is None:
                    continue
                cstart, cend = remap_to_region(cov[0], cov[1], ctx)
                child_entries.append((cstart, cend, ctype))
            child_entries.sort(key=lambda x: (x[0], x[1], x[2]))
            for cstart, cend, ctype in child_entries:
                lines.append(f"{cstart} {cend} {ctype}")
            lines.append("")

    out_path.write_text("\n".join(lines).rstrip() + ("\n" if lines else ""))
    return genes_written


def export_for_mvista(
    species_inputs: List[SpeciesInput],
    reference_species: str,
    mode: str,
    upstream: int,
    downstream: int,
    output_dir: Path,
    export_annotations: bool,
    annotation_features: Set[str],
) -> List[str]:
    output_dir.mkdir(parents=True, exist_ok=True)

    records: List[Tuple[str, str, str, int, int, str]] = []
    log_lines: List[str] = []

    for sp in species_inputs:
        feature = find_gene_feature(sp.gff_path, sp.gene_id)
        region_start, region_end = region_for_feature(feature, mode, upstream, downstream)

        fasta = FastaDB(sp.fasta_path)
        seq = fasta.get_subseq(feature.seqid, region_start, region_end, feature.strand)
        if not seq:
            raise ValueError(
                f"Extracted sequence length is 0: species={sp.species}, region={feature.seqid}:{region_start}-{region_end}"
            )

        header = f"{sp.species}|{sp.gene_id}|{feature.seqid}:{region_start}-{region_end}({feature.strand})"
        records.append((sp.species, header, seq, region_start, region_end, feature.strand))

        sp_fasta = output_dir / f"{sp.species}.mvista.fa"
        sp_fasta.write_text(f">{header}\n{wrap_fasta(seq)}\n")
        log_lines.append(f"Wrote {sp_fasta}")

        if export_annotations:
            output_seqid = safe_seqid(f"{sp.species}|{sp.gene_id}")
            ctx = RegionContext(
                seqid=feature.seqid,
                region_start=region_start,
                region_end=region_end,
                extract_strand=feature.strand,
                output_seqid=output_seqid,
            )
            anno_gff = output_dir / f"{sp.species}.mvista.annotation.gff3"
            anno_bed = output_dir / f"{sp.species}.mvista.annotation.bed"
            anno_txt = output_dir / f"{sp.species}.mvista.annotation.txt"
            gff_records = parse_gff_records(sp.gff_path, feature.seqid)
            region_total, region_child = count_region_records(gff_records, ctx)
            kept = export_remapped_annotations(
                gff_path=sp.gff_path,
                ctx=ctx,
                feature_types=annotation_features,
                gff_out=anno_gff,
                bed_out=anno_bed,
            )
            genes_written = export_mvista_simple_annotation(
                gff_path=sp.gff_path,
                ctx=ctx,
                out_path=anno_txt,
            )
            log_lines.append(
                f"Annotation debug: seqid={feature.seqid}, region_records={region_total}, child_like_records={region_child}"
            )
            log_lines.append(f"Wrote {anno_txt} ({genes_written} genes)")
            log_lines.append(f"Wrote {anno_gff} ({kept} features)")
            log_lines.append(f"Wrote {anno_bed} ({kept} features)")

    # Reference first to simplify upload order.
    records.sort(key=lambda x: (0 if x[0] == reference_species else 1, x[0]))

    multi_fa = output_dir / "mvista_regions.multi.fa"
    with multi_fa.open("w") as out:
        for _sp, header, seq, *_rest in records:
            out.write(f">{header}\n{wrap_fasta(seq)}\n")
    log_lines.append(f"Wrote {multi_fa}")

    manifest = output_dir / "mvista_manifest.tsv"
    with manifest.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["species", "gene_id", "fasta", "annotation_txt", "annotation_gff3", "annotation_bed", "reference"])
        for sp in species_inputs:
            anno_txt_name = f"{sp.species}.mvista.annotation.txt" if export_annotations else ""
            anno_gff_name = f"{sp.species}.mvista.annotation.gff3" if export_annotations else ""
            anno_bed_name = f"{sp.species}.mvista.annotation.bed" if export_annotations else ""
            w.writerow(
                [
                    sp.species,
                    sp.gene_id,
                    f"{sp.species}.mvista.fa",
                    anno_txt_name,
                    anno_gff_name,
                    anno_bed_name,
                    "yes" if sp.species == reference_species else "no",
                ]
            )
    log_lines.append(f"Wrote {manifest}")

    return log_lines


class MVistaGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("mVISTA Ortholog Region Builder")
        self.geometry("1080x640")

        self.species_rows: List[SpeciesInput] = []

        self._build_ui()

    def _build_ui(self) -> None:
        top = ttk.Frame(self, padding=10)
        top.pack(fill="both", expand=True)

        table_frame = ttk.LabelFrame(top, text="Species Inputs", padding=8)
        table_frame.pack(fill="both", expand=True)

        cols = ("species", "fasta", "gff", "gene")
        self.tree = ttk.Treeview(table_frame, columns=cols, show="headings", height=12)
        self.tree.heading("species", text="Species")
        self.tree.heading("fasta", text="Genome FASTA")
        self.tree.heading("gff", text="GFF/GTF")
        self.tree.heading("gene", text="Ortholog Gene ID")
        self.tree.column("species", width=140)
        self.tree.column("fasta", width=320)
        self.tree.column("gff", width=320)
        self.tree.column("gene", width=180)
        self.tree.pack(side="left", fill="both", expand=True)

        yscroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=yscroll.set)
        yscroll.pack(side="right", fill="y")

        form = ttk.Frame(top, padding=(0, 8))
        form.pack(fill="x")

        self.species_var = tk.StringVar()
        self.fasta_var = tk.StringVar()
        self.gff_var = tk.StringVar()
        self.gene_var = tk.StringVar()

        ttk.Label(form, text="Species").grid(row=0, column=0, sticky="w")
        ttk.Entry(form, textvariable=self.species_var, width=18).grid(row=1, column=0, padx=(0, 8))

        ttk.Label(form, text="Genome FASTA").grid(row=0, column=1, sticky="w")
        ttk.Entry(form, textvariable=self.fasta_var, width=36).grid(row=1, column=1, padx=(0, 4))
        ttk.Button(form, text="Browse", command=self._pick_fasta).grid(row=1, column=2, padx=(0, 8))

        ttk.Label(form, text="GFF/GTF").grid(row=0, column=3, sticky="w")
        ttk.Entry(form, textvariable=self.gff_var, width=36).grid(row=1, column=3, padx=(0, 4))
        ttk.Button(form, text="Browse", command=self._pick_gff).grid(row=1, column=4, padx=(0, 8))

        ttk.Label(form, text="Ortholog Gene ID").grid(row=0, column=5, sticky="w")
        ttk.Entry(form, textvariable=self.gene_var, width=20).grid(row=1, column=5, padx=(0, 8))

        ttk.Button(form, text="Add Row", command=self._add_row).grid(row=1, column=6, padx=(0, 6))
        ttk.Button(form, text="Remove Selected", command=self._remove_selected).grid(row=1, column=7)

        settings = ttk.LabelFrame(top, text="Extraction Settings", padding=8)
        settings.pack(fill="x", pady=(8, 0))

        self.mode_var = tk.StringVar(value="Promoter (upstream from TSS)")
        self.up_var = tk.StringVar(value="2000")
        self.down_var = tk.StringVar(value="500")
        self.ref_var = tk.StringVar()
        self.output_var = tk.StringVar(value=str(Path.cwd() / "mvista_output"))
        self.export_anno_var = tk.BooleanVar(value=True)
        self.anno_types_var = tk.StringVar(value="gene,mRNA,transcript,CDS,exon,five_prime_UTR,three_prime_UTR,UTR")

        ttk.Label(settings, text="Mode").grid(row=0, column=0, sticky="w")
        ttk.Combobox(
            settings,
            textvariable=self.mode_var,
            values=["Promoter (upstream from TSS)", "Gene body + flanks"],
            state="readonly",
            width=28,
        ).grid(row=1, column=0, padx=(0, 8))

        ttk.Label(settings, text="Upstream bp").grid(row=0, column=1, sticky="w")
        ttk.Entry(settings, textvariable=self.up_var, width=10).grid(row=1, column=1, padx=(0, 8))

        ttk.Label(settings, text="Downstream bp").grid(row=0, column=2, sticky="w")
        ttk.Entry(settings, textvariable=self.down_var, width=10).grid(row=1, column=2, padx=(0, 8))

        ttk.Label(settings, text="Reference species").grid(row=0, column=3, sticky="w")
        self.ref_combo = ttk.Combobox(settings, textvariable=self.ref_var, state="readonly", width=18)
        self.ref_combo.grid(row=1, column=3, padx=(0, 8))

        ttk.Label(settings, text="Output dir").grid(row=0, column=4, sticky="w")
        ttk.Entry(settings, textvariable=self.output_var, width=42).grid(row=1, column=4, padx=(0, 4))
        ttk.Button(settings, text="Browse", command=self._pick_output).grid(row=1, column=5, padx=(0, 8))

        ttk.Checkbutton(settings, text="Export remapped annotation", variable=self.export_anno_var).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(8, 0)
        )
        ttk.Label(settings, text="Annotation feature types (comma-separated)").grid(row=2, column=2, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Entry(settings, textvariable=self.anno_types_var, width=52).grid(row=2, column=4, columnspan=2, sticky="we", pady=(8, 0))

        ttk.Button(settings, text="Build mVISTA Files", command=self._run).grid(row=1, column=6, rowspan=2, sticky="ns")

        log_box = ttk.LabelFrame(top, text="Log", padding=8)
        log_box.pack(fill="both", expand=True, pady=(8, 0))
        self.log = tk.Text(log_box, height=10)
        self.log.pack(fill="both", expand=True)

    def _pick_fasta(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All files", "*")])
        if p:
            self.fasta_var.set(p)

    def _pick_gff(self) -> None:
        p = filedialog.askopenfilename(filetypes=[("GFF/GTF", "*.gff *.gff3 *.gtf"), ("All files", "*")])
        if p:
            self.gff_var.set(p)

    def _pick_output(self) -> None:
        p = filedialog.askdirectory()
        if p:
            self.output_var.set(p)

    def _append_log(self, msg: str) -> None:
        self.log.insert("end", msg + "\n")
        self.log.see("end")

    def _refresh_reference_options(self) -> None:
        names = [r.species for r in self.species_rows]
        self.ref_combo["values"] = names
        if names and self.ref_var.get() not in names:
            self.ref_var.set(names[0])
        if not names:
            self.ref_var.set("")

    def _add_row(self) -> None:
        species = self.species_var.get().strip()
        fasta = self.fasta_var.get().strip()
        gff = self.gff_var.get().strip()
        gene = self.gene_var.get().strip()

        if not all([species, fasta, gff, gene]):
            messagebox.showerror("Missing fields", "Species, FASTA, GFF/GTF, and Gene ID are required.")
            return

        row = SpeciesInput(species=species, fasta_path=Path(fasta), gff_path=Path(gff), gene_id=gene)
        self.species_rows.append(row)
        self.tree.insert("", "end", values=(species, fasta, gff, gene))
        self._refresh_reference_options()

        self.species_var.set("")
        self.fasta_var.set("")
        self.gff_var.set("")
        self.gene_var.set("")

    def _remove_selected(self) -> None:
        selected = self.tree.selection()
        if not selected:
            return

        idxs = sorted((self.tree.index(item) for item in selected), reverse=True)
        for idx in idxs:
            self.species_rows.pop(idx)
        for item in selected:
            self.tree.delete(item)
        self._refresh_reference_options()

    def _run(self) -> None:
        if len(self.species_rows) < 2:
            messagebox.showerror("Need multiple species", "Add at least 2 species for conservation comparison.")
            return

        ref = self.ref_var.get().strip()
        if not ref:
            messagebox.showerror("Missing reference", "Select a reference species.")
            return

        try:
            upstream = int(self.up_var.get().strip())
            downstream = int(self.down_var.get().strip())
            if upstream < 0 or downstream < 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid length", "Upstream/Downstream must be non-negative integers.")
            return

        out_dir = Path(self.output_var.get().strip())
        mode = self.mode_var.get().strip()
        export_annotations = self.export_anno_var.get()
        annotation_features: Set[str] = set()
        if export_annotations:
            annotation_features = {x.strip().lower() for x in self.anno_types_var.get().split(",") if x.strip()}
            if not annotation_features:
                messagebox.showerror("Invalid annotation types", "Specify at least one annotation feature type.")
                return

        try:
            logs = export_for_mvista(
                species_inputs=self.species_rows,
                reference_species=ref,
                mode=mode,
                upstream=upstream,
                downstream=downstream,
                output_dir=out_dir,
                export_annotations=export_annotations,
                annotation_features=annotation_features,
            )
        except Exception as exc:
            self._append_log(f"ERROR: {exc}")
            messagebox.showerror("Build failed", str(exc))
            return

        self._append_log("Build completed.")
        for line in logs:
            self._append_log(line)
        self._append_log("Tip: Upload reference first in mVISTA, then others (or use multi-FASTA if accepted).")
        if export_annotations:
            self._append_log("Tip: For mVISTA plain-text annotation, use *.mvista.annotation.txt for each extracted sequence.")
        messagebox.showinfo("Done", f"mVISTA input files were written to:\n{out_dir}")


def main() -> None:
    app = MVistaGui()
    app.mainloop()


if __name__ == "__main__":
    main()
