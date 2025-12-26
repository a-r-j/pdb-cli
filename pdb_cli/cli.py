import argparse
import os
import sys

import requests

from . import pdb_summary

try:
    from rich.console import Console
    from rich.table import Table
except ImportError:
    Console = None
    Table = None


def download_pdb(pdb_id, cif=False, outdir="."):
    pdb_id = pdb_id.lower()
    if cif:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        ext = "cif"
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        ext = "pdb"
    out_path = os.path.join(outdir, f"{pdb_id}.{ext}")
    r = requests.get(url)
    if r.status_code == 200:
        with open(out_path, "wb") as f:
            f.write(r.content)
        print(f"Downloaded {pdb_id}.{ext} to {out_path}")
    else:
        print(f"Failed to download {pdb_id}.{ext} (HTTP {r.status_code})")


def download_alphafold(uniprot_id, cif=False, outdir="."):
    """Download AlphaFold prediction from AlphaFold Database"""
    uniprot_id = uniprot_id.upper()
    if cif:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.cif"
        out_path = os.path.join(outdir, f"AF-{uniprot_id}-F1-model_v6.cif")
        ext = "cif"
    else:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
        out_path = os.path.join(outdir, f"AF-{uniprot_id}-F1-model_v6.pdb")
        ext = "pdb"

    r = requests.get(url)
    if r.status_code == 200:
        with open(out_path, "wb") as f:
            f.write(r.content)
        print(
            f"Downloaded AlphaFold prediction for {uniprot_id} ({ext.upper()}) to {out_path}"
        )
    else:
        print(
            f"Failed to download AlphaFold prediction for {uniprot_id} (HTTP {r.status_code})"
        )
        print("Note: AlphaFold predictions may not be available for all UniProt IDs")


def main():
    # Check if the first argument is 'download'
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        parser = argparse.ArgumentParser(
            description="Download a PDB or mmCIF file from RCSB"
        )
        parser.add_argument("download", help=argparse.SUPPRESS)
        parser.add_argument("pdb_id", help="PDB ID to download")
        parser.add_argument(
            "--cif",
            action="store_true",
            help="Download as mmCIF format (default: PDB format)",
        )
        parser.add_argument(
            "--outdir",
            default=".",
            help="Output directory (default: current directory)",
        )
        args = parser.parse_args()
        download_pdb(args.pdb_id, cif=args.cif, outdir=args.outdir)
        return

    # Check if the first argument is 'alphafold'
    if len(sys.argv) > 1 and sys.argv[1] == "alphafold":
        parser = argparse.ArgumentParser(
            description="Download AlphaFold prediction from AlphaFold Database"
        )
        parser.add_argument("alphafold", help=argparse.SUPPRESS)
        parser.add_argument(
            "uniprot_id", help="UniProt ID to download AlphaFold prediction for"
        )
        parser.add_argument(
            "--cif",
            action="store_true",
            help="Download as mmCIF format (default: PDB format)",
        )
        parser.add_argument(
            "--outdir",
            default=".",
            help="Output directory (default: current directory)",
        )
        args = parser.parse_args()
        download_alphafold(args.uniprot_id, cif=args.cif, outdir=args.outdir)
        return

    # Otherwise, treat as summary/analysis
    parser = argparse.ArgumentParser(
        description="Extract summary information from PDB and mmCIF files."
    )
    parser.add_argument("pdb_file", help="Path to the PDB or mmCIF file")
    parser.add_argument(
        "--chains", action="store_true", help="Show number of chains and chain names"
    )
    parser.add_argument(
        "--sequence", action="store_true", help="Show sequence(s) with chain names"
    )
    parser.add_argument(
        "--composition", action="store_true", help="Show residue composition"
    )
    parser.add_argument(
        "--residues", action="store_true", help="Show all residues (chain:code:number)"
    )
    parser.add_argument("--atoms", action="store_true", help="Show number of atoms")
    parser.add_argument(
        "--unique-residues", action="store_true", help="Show number of unique residues"
    )
    parser.add_argument("--bbox", action="store_true", help="Show bounding box")
    parser.add_argument(
        "--ligands", action="store_true", help="Show non-protein ligands"
    )
    parser.add_argument(
        "--ptm", action="store_true", help="Show post-translational modifications"
    )
    parser.add_argument(
        "--missing-residues", action="store_true", help="Show missing residues"
    )
    parser.add_argument(
        "--bfactors", action="store_true", help="Show B-factor statistics"
    )
    parser.add_argument(
        "--binding-sites", action="store_true", help="Show substrate binding sites"
    )
    parser.add_argument(
        "--binding-radius",
        type=float,
        default=5.0,
        help="Radius (Å) for binding site analysis (default: 5.0)",
    )
    parser.add_argument(
        "--interface", action="store_true", help="Show protein-protein interfaces"
    )
    parser.add_argument(
        "--interface-radius",
        type=float,
        default=5.0,
        help="Radius (Å) for interface analysis (default: 5.0)",
    )
    parser.add_argument(
        "--color", action="store_true", help="Enable colored output (requires rich)"
    )

    args = parser.parse_args()

    # Set default behavior: if no flags are provided, show chains and sequence
    if not any(
        [
            args.chains,
            args.sequence,
            args.composition,
            args.residues,
            args.atoms,
            args.unique_residues,
            args.bbox,
            args.ligands,
            args.ptm,
            args.missing_residues,
            args.bfactors,
            args.binding_sites,
            args.interface,
        ]
    ):
        args.chains = True
        args.sequence = True

    if args.color and Console is None:
        print(
            "[ERROR] 'rich' is required for colored output. Please install it.",
            file=sys.stderr,
        )
        sys.exit(1)

    summary = pdb_summary.summarize(
        args.pdb_file,
        chains=args.chains,
        sequence=args.sequence,
        composition=args.composition,
        residues=args.residues,
        atoms=args.atoms,
        unique_residues=args.unique_residues,
        bbox=args.bbox,
        ligands=args.ligands,
        ptm=args.ptm,
        missing_residues=args.missing_residues,
        bfactors=args.bfactors,
        binding_sites=args.binding_sites,
        binding_radius=args.binding_radius,
        interface=args.interface,
        interface_radius=args.interface_radius,
    )

    if args.color:
        console = Console()
        table = Table(
            title=f"[bold blue]Summary for {args.pdb_file}[/bold blue]",
            title_style="bold blue",
            border_style="bright_blue",
            header_style="bold cyan",
            row_styles=["", "dim"],
        )
        table.add_column("Feature", style="bold cyan", no_wrap=True, width=30)
        table.add_column("Value", style="white", no_wrap=False)
        for key, value in summary.items():
            formatted_value = format_value(value, color=True)
            table.add_row(str(key), formatted_value)
        console.print(table)
    else:
        print(f"\n📊 Summary for {args.pdb_file}:")
        print("=" * (len(f"Summary for {args.pdb_file}:") + 15))
        for i, (key, value) in enumerate(summary.items()):
            formatted_value = format_value(value, color=False)
            # Add visual separators for better readability
            if i > 0:
                print("-" * 50)
            print(f"🔹 {key}: {formatted_value}")
        print("=" * (len(f"Summary for {args.pdb_file}:") + 15))


def format_value(value, color=False):
    """Format complex values for better display"""
    if color:
        from rich.console import Group
        from rich.text import Text

    if isinstance(value, dict):
        if color:
            formatted_lines = []
            for k, v in value.items():
                # Create properly styled Text objects
                if "Chain" in k:
                    key_text = Text(k, style="bold yellow")
                    value_text = Text(str(v), style="green")
                elif "Count" in k or "Number" in k:
                    key_text = Text(k, style="bold magenta")
                    value_text = Text(str(v), style="bright_cyan")
                elif "Distance" in str(v):
                    key_text = Text(k, style="bold red")
                    value_text = Text(str(v), style="bright_green")
                elif "B-factor" in k:
                    key_text = Text(k, style="bold orange3")
                    value_text = Text(str(v), style="bright_yellow")
                elif "ligand" in k.lower() or "binding" in k.lower():
                    key_text = Text(k, style="bold purple")
                    value_text = Text(str(v), style="bright_magenta")
                elif "binding" in k.lower() or "ligand" in k.lower():
                    key_text = Text(k, style="bold purple")
                    value_text = Text(str(v), style="bright_magenta")
                elif "Sequence" in k:
                    key_text = Text(k, style="bold blue")
                    # For sequences, just return the raw string to let Rich handle wrapping
                    if isinstance(v, dict):
                        sequence_text = ""
                        for chain_key, sequence_value in v.items():
                            sequence_text += f"{chain_key}: {sequence_value}\n"
                        value_text = sequence_text.rstrip()
                    else:
                        value_text = str(v)
                elif "Residues" in k and isinstance(v, list):
                    key_text = Text(k, style="bold green")
                    # For residues list, format as a readable list
                    if len(v) > 20:
                        # Show first 10 and last 10 if list is long
                        first_part = v[:10]
                        last_part = v[-10:]
                        residue_text = (
                            "\n".join(first_part)
                            + f"\n... ({len(v) - 20} more residues) ...\n"
                            + "\n".join(last_part)
                        )
                    else:
                        residue_text = "\n".join(v)
                    value_text = residue_text
                else:
                    key_text = Text(k, style="bold white")
                    value_text = Text(str(v), style="green")

                # Combine key and value with separator
                line = key_text + Text(": ", style="white") + value_text
                formatted_lines.append(line)

            return Group(*formatted_lines)
        else:
            formatted = []
            for k, v in value.items():
                if "Sequence" in k and isinstance(v, dict):
                    formatted.append(f"  {k}:")
                    for chain_key, sequence_value in v.items():
                        sequence_text = str(sequence_value)
                        if len(sequence_text) > 80:
                            # Split long sequences into chunks
                            chunks = [
                                sequence_text[i : i + 80]
                                for i in range(0, len(sequence_text), 80)
                            ]
                            formatted.append(f"    {chain_key}: {chunks[0]}")
                            for chunk in chunks[1:]:
                                formatted.append(f"      {chunk}")
                        else:
                            formatted.append(f"    {chain_key}: {sequence_text}")
                elif "Residues" in k and isinstance(v, list):
                    formatted.append(f"  {k}:")
                    if len(v) > 20:
                        # Show first 10 and last 10 if list is long
                        first_part = v[:10]
                        last_part = v[-10:]
                        for item in first_part:
                            formatted.append(f"    {item}")
                        formatted.append(f"    ... ({len(v) - 20} more residues) ...")
                        for item in last_part:
                            formatted.append(f"    {item}")
                    else:
                        for item in v:
                            formatted.append(f"    {item}")
                else:
                    formatted.append(f"  {k}: {v}")
            return "\n" + "\n".join(formatted)
    elif isinstance(value, list):
        if color:
            items = [Text(str(item), style="bright_cyan") for item in value]
            result = Text()
            for i, item in enumerate(items):
                if i > 0:
                    result += Text("\n", style="white")
                result += item
            return result
        else:
            return "\n".join(str(item) for item in value)
    elif isinstance(value, str) and "No" in value:
        if color:
            return Text(value, style="dim italic")
        else:
            return value
    elif isinstance(value, str) and ("Å" in value or "distance" in value.lower()):
        if color:
            return Text(value, style="bright_green")
        else:
            return value
    elif isinstance(value, str) and ("Ligand:" in value or "binding" in value.lower()):
        if color:
            # Handle binding site output with custom formatting
            from rich.console import Group
            from rich.text import Text

            lines = value.split("\n")
            formatted_lines = []

            for line in lines:
                if line.strip().startswith("Ligand:"):
                    # Ligand names in bold yellow
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        formatted_lines.append(
                            Text("  Ligand:", style="bold white")
                            + Text(":", style="white")
                            + Text(parts[1].strip(), style="bold yellow")
                        )
                elif "Number of interacting residues:" in line:
                    # Count in bright cyan
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        formatted_lines.append(
                            Text(
                                "    Number of interacting residues:",
                                style="bold white",
                            )
                            + Text(":", style="white")
                            + Text(parts[1].strip(), style="bright_cyan")
                        )
                elif line.strip().startswith("Residues:"):
                    # Section header
                    formatted_lines.append(Text(line, style="bold white"))
                elif line.strip().startswith("      ") and ":" in line:
                    # Residue entries
                    parts = line.split(":", 2)
                    if len(parts) >= 3:
                        residue_part = parts[0].strip() + ":" + parts[1].strip()
                        distance_part = ":" + parts[2]
                        formatted_lines.append(
                            Text("      ", style="white")
                            + Text(residue_part, style="green")
                            + Text(distance_part, style="bright_green")
                        )
                else:
                    # Regular lines
                    formatted_lines.append(Text(line, style="white"))

            return Group(*formatted_lines)
        else:
            return value
    elif isinstance(value, str) and (
        "Interface" in value or "interface" in value.lower()
    ):
        if color:
            # Handle interface output with custom formatting
            from rich.console import Group
            from rich.text import Text

            lines = value.split("\n")
            formatted_lines = []

            for line in lines:
                if line.strip().startswith("Interface ") and ":" not in line:
                    # Interface names in bold cyan
                    parts = line.split(" ", 1)
                    if len(parts) == 2:
                        formatted_lines.append(
                            Text("  ", style="white")
                            + Text("Interface ", style="bold white")
                            + Text(parts[1], style="bold cyan")
                        )
                elif "Number of interface residues:" in line:
                    # Count in bright cyan
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        formatted_lines.append(
                            Text(
                                "    Number of interface residues:",
                                style="bold white",
                            )
                            + Text(":", style="white")
                            + Text(parts[1].strip(), style="bright_cyan")
                        )
                elif line.strip().startswith("Residues:"):
                    # Section header
                    formatted_lines.append(Text(line, style="bold white"))
                elif line.strip().startswith("      ") and ":" in line:
                    # Residue entries
                    parts = line.split(":", 2)
                    if len(parts) >= 3:
                        residue_part = parts[0].strip() + ":" + parts[1].strip()
                        distance_part = ":" + parts[2]
                        formatted_lines.append(
                            Text("      ", style="white")
                            + Text(residue_part, style="green")
                            + Text(distance_part, style="bright_green")
                        )
                else:
                    # Regular lines
                    formatted_lines.append(Text(line, style="white"))

            return Group(*formatted_lines)
        else:
            return value
    else:
        if color:
            return Text(str(value), style="white")
        else:
            return str(value)
