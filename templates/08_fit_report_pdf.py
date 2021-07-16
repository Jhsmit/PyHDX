"""Generate a pdf output with all peptide fits. Requires pdflatex"""
from pyhdx.output import Output, Report
from pyhdx.fileIO import load_fitresult
from pathlib import Path

current_dir = Path().cwd()
fit_result = load_fitresult(current_dir / 'output' / 'SecB_fit')

output = Output(fit_result)

report = Report(output)
report.add_peptide_figures()
report.generate_pdf(current_dir / 'output' / 'SecB_fit_report')