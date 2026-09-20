from datetime import datetime

def extract_date_from_header(header):
    parts = header.split('|')
    if len(parts) < 2:
        return None
    
    date_str = parts[-1].strip()
    
    try:
        if len(date_str) == 10:
            d = datetime.strptime(date_str, "%Y-%m-%d")
            return d.year + (d.timetuple().tm_yday - 1) / 365.0
        elif len(date_str) == 7:
            year, month = map(int, date_str.split('-'))
            return year + ((month - 1) / 12.0) + (15 / 365.0)
        elif len(date_str) == 4:
            return float(date_str) + 0.5
    except:
        pass
    return None

input_fasta = "../resources/colombia_data/aligned/align_DENV1_BEAST.fasta"
output_txt = "../resources/BEAST/DENV1_dates_beauti.txt"

with open(input_fasta, 'r') as f_in, open(output_txt, 'w') as f_out:
    for line in f_in:
        if line.startswith('>'):
            header = line[1:].strip()
            date_decimal = extract_date_from_header(header)
            if date_decimal is not None:
                f_out.write(f"{header}\t{date_decimal:.4f}\n")

print(f"Archivo generado: {output_txt}")