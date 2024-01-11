"""Load two HDX-MS datasets and guesses and perform fitting in batch with a second regualizer and mock alignment"""

raise DeprecationWarning("Outdated example")

from pathlib import Path

from pyhdx.batch_processing import StateParser
from pyhdx.fitting import fit_gibbs_global_batch_aligned
from pyhdx.fileIO import csv_to_protein
import yaml


mock_alignment = {
    "dimer": "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------",
    "apo": "MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA",
}

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)


test_data_dir = current_dir.parent / "tests" / "test_data"
guess = csv_to_protein(test_data_dir / "output" / "ecSecB_guess.csv")

input_dir = test_data_dir / "input"
parser = StateParser(hdx_spec, input_dir)


hdx_set = hdxm = parser.load_hdxmset()
gibbs_guess = hdx_set[0].guess_deltaG(guess["rate"])

hdx_set.add_alignment(list(mock_alignment.values()))
result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)

# Human readable output
result.to_file(output_dir / "Batch_aligned_fit_result.txt", fmt="pprint")

# Machine readable output
result.to_file(output_dir / "Batch_aligned_fit_result.csv", fmt="csv")
