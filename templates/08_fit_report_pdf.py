"""Generate a pdf output with all peptide fits. Requires pdflatex"""
from pyhdx.output import FitReport
from pyhdx.fileIO import load_fitresult
from pathlib import Path
from concurrent import futures

current_dir = Path().cwd()
fit_result = load_fitresult(current_dir / 'output' / 'SecB_tetramer_dimer_batch')

tmp_dir = Path(__file__).parent / 'temp'
tmp_dir.mkdir(exist_ok=True)

if __name__ == '__main__':

    report = FitReport(fit_result, temp_dir=tmp_dir)
    report.add_peptide_uptake_curves()

    executor = futures.ProcessPoolExecutor(max_workers=10)

    report.generate_figures(executor=executor)
    report.generate_latex()
    report.generate_pdf(current_dir / 'pdftest123')