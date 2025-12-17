from hdxms_datasets.process import merge_peptide_tables, compute_uptake_metrics
from hdxms_datasets.utils import get_peptides_by_type
from hdxms_datasets.models import State, DeuterationType, Peptides
import narwhals as nw
import pandas as pd
import warnings

from pyhdx.process import correct_d_uptake


def state_kwargs(state: State):
    return {
        "name": state.name,
        "sequence": state.protein_state.sequence,
        "n_term": state.protein_state.n_term,
        "c_term": state.protein_state.c_term,
    }


def peptides_kwargs(peptides: Peptides):
    return {
        "temperature": peptides.temperature,
        "pH": peptides.pH,
    }


def adapt_for_pyhdx(df: nw.DataFrame) -> nw.DataFrame:
    """adapt open hdx dataframes to match pyhdx expectations"""

    rename = {
        "frac_fd_control": "rfu",
        "frac_fd_control_sd": "rfu_sd",
    }

    df = (
        df.rename(rename).with_columns((nw.col("end") + 1).alias("stop")).drop_nulls(subset=["rfu"])
    )
    return df


def load_peptides(peptides):
    deuteration_types = [
        DeuterationType.partially_deuterated,
        DeuterationType.fully_deuterated,
        DeuterationType.non_deuterated,
    ]

    output = {}
    for dtype in deuteration_types:
        p = get_peptides_by_type(peptides, dtype)
        if p is None:
            continue
        output[dtype.name] = p.load()

    return output


def parse_dataset_states(
    states: list[State], drop_first: int, d_percentage: float | None = None
) -> list[tuple[pd.DataFrame, dict]]:
    """Parse an HDXDataSet into a list of tuples of (peptides, metadata) for pyhdx"""

    output = []
    loaded_peptides = [load_peptides(state.peptides) for state in states]

    def find_fd_control() -> tuple[int, pd.DataFrame]:
        for i, lp in enumerate(loaded_peptides):
            if "fully_deuterated" in lp:
                return i, lp["fully_deuterated"]
        raise ValueError("No FD control found")

    if all("fully_deuterated" in lp for lp in loaded_peptides):
        pass
    else:
        idx, fd_control = find_fd_control()
        warnings.warn(
            f"Not all states have FD control, using FD control from state index {idx}: {states[idx].name}"
        )
        # add the fd_control to all loaded peptides that lack it
        for lp in loaded_peptides:
            if "fully_deuterated" not in lp:
                lp["fully_deuterated"] = fd_control

    for state, peptides in zip(states, loaded_peptides):
        pd_peptides = get_peptides_by_type(state.peptides, DeuterationType.partially_deuterated)
        assert pd_peptides is not None  # this never happens due to previous checks

        d_percentage = d_percentage or pd_peptides.d_percentage
        assert d_percentage is not None, (
            "Deuterium percentage must be specified either in the dataset or as a kwarg"
        )
        merged = merge_peptide_tables(**peptides)  # type: ignore
        computed = compute_uptake_metrics(merged)
        adapted = adapt_for_pyhdx(computed).to_pandas()
        peptides_corrected = correct_d_uptake(
            adapted,
            drop_first=drop_first,
            d_percentage=d_percentage,
        )

        metadata = {
            **state_kwargs(state),
            **peptides_kwargs(pd_peptides),
        }

        output.append((peptides_corrected, metadata))

    return output
