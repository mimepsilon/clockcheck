"""Command line: `clockcheck audit ...`, `clockcheck demo`, `clockcheck list`."""
from __future__ import annotations
import argparse
import sys

from .clocks import BUILTIN, LinearClock
from .data import load_csv, load_geo, synthetic_fixture
from .report import build_audit_card, card_to_json, card_to_markdown


def _run_audit(args) -> int:
    clock = LinearClock.load(args.clock)
    if args.geo:
        ds = load_geo(args.geo, age_col=args.age_col)
    elif args.betas and args.meta:
        ds = load_csv(args.betas, args.meta, age_col=args.age_col)
    else:
        print("Provide --geo GSE_ID or both --betas and --meta.", file=sys.stderr)
        return 2
    compare = [LinearClock.load(n) for n in (args.compare or []) if n != args.clock]
    card = build_audit_card(clock, ds, compare_clocks=compare or None)
    out = card_to_json(card) if args.format == "json" else card_to_markdown(card)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(out)
        print(f"wrote {args.out}")
    else:
        print(out)
    return 0


def _run_demo(_args) -> int:
    clocks = [LinearClock.load("Horvath1"), LinearClock.load("Hannum")]
    ds = synthetic_fixture(clocks, n_samples=60, seed=0)
    card = build_audit_card(clocks[0], ds, compare_clocks=clocks[1:])
    print(card_to_markdown(card))
    return 0


def _run_list(_args) -> int:
    for name, spec in BUILTIN.items():
        print(f"{name:12s} unit={spec.unit:6s} {spec.citation}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="clockcheck",
                                description="Honest audits for epigenetic aging clocks.")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("audit", help="audit a clock on a dataset")
    a.add_argument("--clock", required=True, choices=list(BUILTIN))
    a.add_argument("--geo", help="GEO series id, e.g. GSE40279 (needs biolearn + network)")
    a.add_argument("--betas", help="CpG x sample beta matrix CSV")
    a.add_argument("--meta", help="sample metadata CSV (must include the age column)")
    a.add_argument("--age-col", default="age")
    a.add_argument("--compare", nargs="*", help="other clocks for cross-clock agreement")
    a.add_argument("--format", choices=["md", "json"], default="md")
    a.add_argument("--out", help="write to file instead of stdout")
    a.set_defaults(func=_run_audit)

    sub.add_parser("demo", help="run on the labeled synthetic fixture (offline)"
                   ).set_defaults(func=_run_demo)
    sub.add_parser("list", help="list built-in clocks").set_defaults(func=_run_list)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
