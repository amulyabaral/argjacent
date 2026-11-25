#!/usr/bin/env python3
"""Compacted gene-centered local GFA v1.0 exporter for ARGjacent.

Builds a lossless, blunt-end (`0M` overlap) directed acyclic graph per gene:

    [upstream left suffix tree] ---> [GENE (5'->3')] ---> [downstream right prefix tree]

Properties:
  1. Lossless bijection: every candidate context in `flanks.json` corresponds to a unique
     source->GENE (left) or GENE->sink (right) walk, emitted as a GFA `P` path record.
  2. Explicit certainty split: any branch segment crossed by `context_len` (`certain_bp`) or
     `end_len` (`certain_end_bp`) is split at the exact certified boundary so every `S` segment
     carries an unambiguous tier tag and Bandage color:
       - `TP:Z:gene`      (`CL:Z:#d95f02` orange) : focal gene (5'->3')
       - `TP:Z:certain`   (`CL:Z:#1b9e77` teal)   : certified flank (`context_len`)
       - `TP:Z:end_only`  (`CL:Z:#66a61e` green)  : end-anchored extension (`end_len` > `context_len`)
       - `TP:Z:candidate` (`CL:Z:#b3b3b3` grey)   : unverified candidate extension beyond certainty
  3. Local segment depths (`DP:f:`): computed from `depth_bins` (100-bp binned exact depth from
     `verify.py`) over each segment's exact outward interval `[d0, d1)`, falling back to `depth`
     for legacy 0.3/0.4 JSON files.
  4. Optional two-sided `P` lines (`<gene>|pair|p1|<tier>`) and cohort repeat-resolved `P` lines.
"""
import argparse
import json
import os
import statistics

DEPTH_BIN = 100


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, encoding="latin-1") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf).upper()))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            recs.append((name, "".join(buf).upper()))
    return recs


def _sanitize_id(s):
    """GFA segment/path IDs cannot contain whitespace or GFA path delimiters (`,`, `;`, `+`, `-`)."""
    out = []
    for ch in s:
        out.append("_" if ch in " \t,;+-" else ch)
    return "".join(out)


def _segment_depth(cands, d0, d1, fallback_dp):
    """Median local depth across outward coordinates `[d0, d1)` using `depth_bins` if available."""
    vals = []
    b0 = d0 // DEPTH_BIN
    b1 = max(b0 + 1, (d1 + DEPTH_BIN - 1) // DEPTH_BIN)
    for p in cands:
        bins = p.get("depth_bins")
        if bins and b0 < len(bins):
            sl = bins[b0:min(len(bins), b1)]
            if sl:
                vals.append(float(statistics.median(sl)))
    if vals:
        return round(max(vals), 1)
    return round(fallback_dp, 1)


def _build_outward_segments(paths, side, gid_prefix):
    """Build compacted outward segments from the gene boundary for `side` ('left' or 'right')."""
    if not paths:
        return [], [], {}, []

    ch_map = [{}]
    node_cands = [[]]
    node_dist = [0]
    edge_tier = [0]
    edge_depth = [0.0]
    edge_main = [False]
    end_here = [set()]

    cand_end_node = {}
    for ci, p in enumerate(paths):
        flk = p["flank"]
        if not flk:
            end_here[0].add(ci)
            cand_end_node[ci] = 0
            continue
        seq_out = flk if side == "right" else flk[::-1]
        cert = int(p.get("context_len", 0))
        endc = max(cert, int(p.get("end_len", cert)))
        dp = float(p.get("depth", 0.0))
        is_m = bool(p.get("main", False))
        u = 0
        node_cands[0].append(ci)
        for d, base in enumerate(seq_out):
            tier = 2 if d < cert else (1 if d < endc else 0)
            nxt = ch_map[u].get(base)
            if nxt is None:
                nxt = len(ch_map)
                ch_map[u][base] = nxt
                ch_map.append({})
                node_cands.append([])
                node_dist.append(d + 1)
                edge_tier.append(tier)
                edge_depth.append(dp)
                edge_main.append(is_m)
                end_here.append(set())
            else:
                if tier > edge_tier[nxt]:
                    edge_tier[nxt] = tier
                if dp > edge_depth[nxt]:
                    edge_depth[nxt] = dp
                if is_m:
                    edge_main[nxt] = True
            u = nxt
            node_cands[u].append(ci)
        end_here[u].add(ci)
        cand_end_node[ci] = u

    colors = {2: ("certain", "#1b9e77"), 1: ("end_only", "#66a61e"), 0: ("candidate", "#b3b3b3")}
    seg_idx = 0
    raw_segs = []
    node_to_out_segs = {0: []}

    stack = [(0, b, v) for b, v in sorted(ch_map[0].items(), reverse=True)]
    while stack:
        u0, b0, v = stack.pop()
        chars = [b0]
        tier = edge_tier[v]
        dp = edge_depth[v]
        is_m = edge_main[v]
        cset = list(node_cands[v])
        cur = v
        while len(ch_map[cur]) == 1 and not end_here[cur]:
            (b_next, nxt), = ch_map[cur].items()
            if edge_tier[nxt] != tier:
                break
            chars.append(b_next)
            if edge_depth[nxt] > dp:
                dp = edge_depth[nxt]
            if edge_main[nxt]:
                is_m = True
            cur = nxt
        seg_idx += 1
        sid = "%s|%s%d" % (gid_prefix, "L" if side == "left" else "R", seg_idx)
        out_str = "".join(chars)
        seq_5to3 = out_str if side == "right" else out_str[::-1]
        tp, cl = colors[tier]
        d0 = node_dist[u0]
        d1 = node_dist[cur]
        seg_dp = _segment_depth([paths[c] for c in cset], d0, d1, dp)
        s_idx = len(raw_segs)
        raw_segs.append({
            "id": sid, "seq": seq_5to3, "ln": len(seq_5to3), "dp": seg_dp,
            "tp": tp, "cl": cl, "main": is_m, "u0": u0, "u1": cur,
            "d0": d0, "d1": d1,
            "cands": [paths[c]["id"] for c in cset],
        })
        node_to_out_segs.setdefault(u0, []).append(s_idx)
        node_to_out_segs.setdefault(cur, [])
        for b_next, nxt in sorted(ch_map[cur].items(), reverse=True):
            stack.append((cur, b_next, nxt))

    u1_to_seg = {s["u1"]: i for i, s in enumerate(raw_segs)}
    root_seg_ids = [raw_segs[i]["id"] for i in node_to_out_segs.get(0, [])]

    links = []
    for s in raw_segs:
        for child_i in node_to_out_segs.get(s["u1"], []):
            cs = raw_segs[child_i]
            if side == "right":
                links.append((s["id"], "+", cs["id"], "+", "0M"))
            else:
                links.append((cs["id"], "+", s["id"], "+", "0M"))

    cand_seg_lists = {}
    for ci, p in enumerate(paths):
        cid = p["id"]
        u = cand_end_node[ci]
        chain = []
        while u != 0:
            si = u1_to_seg[u]
            chain.append(raw_segs[si]["id"])
            u = raw_segs[si]["u0"]
        cand_seg_lists[cid] = chain if side == "left" else chain[::-1]

    return raw_segs, links, cand_seg_lists, root_seg_ids


def build_gene_gfa(gene_name, gene_seq, gene_entry, pairs=None, cohort_exts=None,
                   gid_override=None, read_len=150.0):
    """Build GFA v1.0 lines (`S`, `L`, `P`) for one gene.

    Returns `(lines, seg_seq, l_cpaths, r_cpaths)`.
    """
    gid = gid_override if gid_override is not None else _sanitize_id(gene_name)
    rlen = max(1.0, float(read_len) if read_len else 150.0)
    left_paths = gene_entry.get("left", {}).get("paths", [])
    right_paths = gene_entry.get("right", {}).get("paths", [])
    g_deps = [float(p.get("gene_depth", 0.0)) for p in left_paths + right_paths if p.get("gene_depth", 0.0) > 0]
    if g_deps:
        g_dp = round(max(g_deps), 1)
    else:
        g_dp = round(max(
            float(gene_entry.get("left", {}).get("read_depth", 0.0)),
            float(gene_entry.get("right", {}).get("read_depth", 0.0)),
        ), 1)

    l_segs, l_links, l_cpaths, l_roots = _build_outward_segments(left_paths, "left", gid)
    r_segs, r_links, r_cpaths, r_roots = _build_outward_segments(right_paths, "right", gid)

    g_seg_id = "%s|GENE" % gid
    lines = []
    seg_seq = {g_seg_id: gene_seq}

    for s in l_segs:
        seg_seq[s["id"]] = s["seq"]
        rc = int(round(s["dp"] * s["ln"] / rlen))
        lines.append("S\t%s\t%s\tLN:i:%d\tDP:f:%.1f\tRC:i:%d\tTP:Z:%s\tSD:Z:left\tCL:Z:%s" % (
            s["id"], s["seq"], s["ln"], s["dp"], rc, s["tp"], s["cl"]))
    rc_g = int(round(g_dp * len(gene_seq) / rlen))
    lines.append("S\t%s\t%s\tLN:i:%d\tDP:f:%.1f\tRC:i:%d\tTP:Z:gene\tSD:Z:gene\tCL:Z:#d95f02" % (
        g_seg_id, gene_seq, len(gene_seq), g_dp, rc_g))
    for s in r_segs:
        seg_seq[s["id"]] = s["seq"]
        rc = int(round(s["dp"] * s["ln"] / rlen))
        lines.append("S\t%s\t%s\tLN:i:%d\tDP:f:%.1f\tRC:i:%d\tTP:Z:%s\tSD:Z:right\tCL:Z:%s" % (
            s["id"], s["seq"], s["ln"], s["dp"], rc, s["tp"], s["cl"]))

    all_links = list(l_links)
    for lr in l_roots:
        all_links.append((lr, "+", g_seg_id, "+", "0M"))
    for rr in r_roots:
        all_links.append((g_seg_id, "+", rr, "+", "0M"))
    all_links.extend(r_links)
    for u, uo, v, vo, cig in all_links:
        lines.append("L\t%s\t%s\t%s\t%s\t%s" % (u, uo, v, vo, cig))

    for p in left_paths:
        cid = p["id"]
        slist = l_cpaths[cid] + [g_seg_id]
        lines.append("P\t%s|left|%s\t%s\t*" % (gid, cid, ",".join(x + "+" for x in slist)))
    for p in right_paths:
        cid = p["id"]
        slist = [g_seg_id] + r_cpaths[cid]
        lines.append("P\t%s|right|%s\t%s\t*" % (gid, cid, ",".join(x + "+" for x in slist)))

    if pairs:
        for pr in pairs:
            lid, rid = pr["left_id"], pr["right_id"]
            if lid in l_cpaths and rid in r_cpaths:
                slist = l_cpaths[lid] + [g_seg_id] + r_cpaths[rid]
                lines.append("P\t%s|pair|%s|%s\t%s\t*" % (gid, pr["id"], pr["tier"], ",".join(x + "+" for x in slist)))

    if cohort_exts:
        for sd_name, exts in sorted(cohort_exts.items()):
            cpaths = l_cpaths if sd_name == "left" else r_cpaths
            for ex in exts:
                cid = ex.get("cand_id")
                if cid in cpaths:
                    slist = (cpaths[cid] + [g_seg_id]) if sd_name == "left" else ([g_seg_id] + cpaths[cid])
                    lines.append("P\t%s|cohort_ext|%s|%s|%dbp\t%s\t*" % (
                        gid, sd_name, cid, int(ex.get("cohort_certified_bp", 0)), ",".join(x + "+" for x in slist)))

    return lines, seg_seq, l_cpaths, r_cpaths


def verify_roundtrip(gene_entry, seg_seq, l_cpaths, r_cpaths):
    """Assert 100 % byte-for-byte sequence reconstruction of every candidate path."""
    for p in gene_entry.get("left", {}).get("paths", []):
        recon = "".join(seg_seq[s] for s in l_cpaths[p["id"]])
        if recon != p["flank"]:
            raise AssertionError("Left GFA roundtrip mismatch on %s (%d vs %d bp)" % (p["id"], len(recon), len(p["flank"])))
    for p in gene_entry.get("right", {}).get("paths", []):
        recon = "".join(seg_seq[s] for s in r_cpaths[p["id"]])
        if recon != p["flank"]:
            raise AssertionError("Right GFA roundtrip mismatch on %s (%d vs %d bp)" % (p["id"], len(recon), len(p["flank"])))


def write_gfa(genes_dict, doc_genes, out_gfa, pairs_by_gene=None, cohort_exts_by_gene=None,
              verify=True, read_len=150.0):
    """Write a combined GFA v1.0 file from in-memory gene sequences and `flanks.json['genes']`."""
    out_lines = ["H\tVN:Z:1.0"]
    n_genes = n_segs = n_links = n_paths = 0
    raw_bp = comp_bp = 0
    used_gids = set()
    for gname, gentry in doc_genes.items():
        gseq = genes_dict.get(gname)
        if not gseq:
            continue
        base_gid = _sanitize_id(gname)
        gid = base_gid
        suffix = 2
        while gid in used_gids:
            gid = "%s__%d" % (base_gid, suffix)
            suffix += 1
        used_gids.add(gid)
        n_genes += 1
        prs = (pairs_by_gene or {}).get(gname)
        if prs is None:
            prs = gentry.get("pairs")
        cext = (cohort_exts_by_gene or {}).get(gname)
        lines, seg_seq, lc, rc = build_gene_gfa(
            gname, gseq, gentry, pairs=prs, cohort_exts=cext,
            gid_override=gid, read_len=read_len,
        )
        if verify:
            verify_roundtrip(gentry, seg_seq, lc, rc)
        for ln in lines:
            if ln.startswith("S\t"):
                n_segs += 1
                comp_bp += len(ln.split("\t")[2])
            elif ln.startswith("L\t"):
                n_links += 1
            elif ln.startswith("P\t"):
                n_paths += 1
        raw_bp += len(gseq) + sum(len(p["flank"]) for sd in ("left", "right") for p in gentry.get(sd, {}).get("paths", []))
        out_lines.extend(lines)
    parent = os.path.dirname(os.path.abspath(out_gfa))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(out_gfa, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    return {
        "genes": n_genes, "segments": n_segs, "links": n_links, "paths": n_paths,
        "raw_bp": raw_bp, "compacted_bp": comp_bp,
        "compression_ratio": round(raw_bp / max(1, comp_bp), 2),
    }


def export_sample_gfa(genes_fa, flanks_json, out_gfa, pairs_by_gene=None, cohort_exts_by_gene=None, verify=True):
    genes = dict(read_fasta(genes_fa))
    doc = json.load(open(flanks_json))
    rlen = doc.get("parameters", {}).get("read_len") or 150.0
    return write_gfa(genes, doc.get("genes", {}), out_gfa, pairs_by_gene=pairs_by_gene,
                     cohort_exts_by_gene=cohort_exts_by_gene, verify=verify, read_len=rlen)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Export compacted gene-centered local GFA v1.0 from ARGjacent outputs.")
    ap.add_argument("-g", "--genes", required=True, help="genes.fa input")
    ap.add_argument("-j", "--json", required=True, help="PREFIX.flanks.json")
    ap.add_argument("-o", "--out", required=True, help="output .gfa file")
    args = ap.parse_args(argv)
    st = export_sample_gfa(args.genes, args.json, args.out, verify=True)
    print("Wrote %s: %d genes, %d segments, %d links, %d paths (%d -> %d bp, %.2fx compaction)" % (
        args.out, st["genes"], st["segments"], st["links"], st["paths"], st["raw_bp"], st["compacted_bp"], st["compression_ratio"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
