import os

DEFAULT_CATEGORIES = ('CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA')

FASTA_EXTS = {'.fasta', '.fna', '.fa'}
GFF_EXTS = {'.gff', '.gff3'}
GENBANK_EXTS = {'.gbk', '.gb', '.gbff'}


def load_reference(reference_path, companion_path=None, categories=DEFAULT_CATEGORIES):
    '''
    Load a reference genome and its feature annotations.

    Accepts:
        - GenBank  (.gbk / .gb / .gbff)   — sequence and features included.
        - FASTA    (.fasta / .fna / .fa)  — requires companion GFF
        - GFF      (.gff / .gff3)         - requires companion FASTA or embedded FASTA sequence
    '''

    ext = os.path.splitext(reference_path)[1].lower()

    if ext in GENBANK_EXTS:
        return load_genbank(reference_path, categories)

    if ext in FASTA_EXTS:
        if not companion_path:
            raise ValueError(
                f'FASTA reference {reference_path} requires a companion --gff for annotations.'
            )
        seq = load_fasta_file(reference_path)
        gff = load_gff_file(companion_path, categories)
        return seq, gff

    if ext in GFF_EXTS:
        seq, gff = load_gff_with_optional_fasta(reference_path, categories)
        if not seq:
            if not companion_path:
                raise ValueError(
                    f'GFF reference {reference_path} has no embedded FASTA.'
                    'Provide a companion --fasta.'
                )
            seq = load_fasta_file(companion_path)
        return seq, gff

    raise ValueError(f'Unsupported reference extension: {ext} (path: {reference_path})')


def load_sequences(reference_path, companion_path=None):
    '''Return only the contig to sequence map (no annotations needed).'''
    ext = os.path.splitext(reference_path)[1].lower()
    if ext in GENBANK_EXTS:
        return load_genbank(reference_path, DEFAULT_CATEGORIES)[0]
    if ext in FASTA_EXTS:
        return load_fasta_file(reference_path)
    if ext in GFF_EXTS:
        seq, _ = load_gff_with_optional_fasta(reference_path, DEFAULT_CATEGORIES)
        if not seq:
            if not companion_path:
                raise ValueError(
                    f'GFF reference {reference_path} has no embedded FASTA.'
                    'Provide a companion FASTA.')
            seq = load_fasta_file(companion_path)
        return seq
    raise ValueError(f'Unsupported reference extension: {ext} (path: {reference_path})')


def load_fasta_file(path):
    with open(path, 'r') as f:
        return parse_fasta(f)


def parse_fasta(lines):
    seq = {}
    current = None
    buf = []
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current is not None:
                seq[current] = ''.join(buf).upper()
            current = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if current is not None:
        seq[current] = ''.join(buf).upper()
    return seq


def load_gff_file(path, categories):
    gff = {}

    with open(path, 'r') as f:
        for line in f:
            if line.startswith('##FASTA'):
                break
            ingest_gff_line(line, gff, categories)
    return gff


def load_gff_with_optional_fasta(path, categories):
    gff = {}
    seq = {}

    with open(path, 'r') as f:
        in_fasta = False
        fasta_lines = []
        for line in f:
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line.startswith('##FASTA'):
                in_fasta = True
                continue
            ingest_gff_line(line, gff, categories)
        if fasta_lines:
            seq = parse_fasta(fasta_lines)
    return seq, gff


def ingest_gff_line(line, gff, categories):
    if not line or line.startswith('#'):
        return
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 8:
        return
    contig, _source, gene_type, start_s, end_s, _score, strand = parts[:7]
    if gene_type not in categories:
        return
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        return
    name = f'{contig}_{start}_{end}'
    bucket = gff.setdefault(contig, {})
    feature = [name, start, end, strand, gene_type]
    for pos in range(start, end + 1):
        bucket.setdefault(pos, []).append(feature)


def load_genbank(path, categories):
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ImportError(f'{exc}\n'
            'Reading GenBank references requires biopython. Install it with `pip install Bio`.')

    seq = {}
    gff = {}
    allowed = {c.lower() for c in categories}

    for record in SeqIO.parse(path, 'genbank'):
        contig = record.id
        seq[contig] = str(record.seq).upper()
        bucket = gff.setdefault(contig, {})
        for feat in record.features:
            ftype = feat.type
            if ftype.lower() not in allowed:
                continue

            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            strand = '+' if feat.location.strand in (1, None) else '-'
            name = f'{contig}_{start}_{end}'

            # Normalize category names to match defined categories
            canonical = {c.lower(): c for c in categories}[ftype.lower()]
            feature = [name, start, end, strand, canonical]
            for pos in range(start, end + 1):
                bucket.setdefault(pos, []).append(feature)
    return seq, gff
