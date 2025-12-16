from hdxms_datasets.utils import get_peptides_by_type
from hdxms_datasets.models import HDXDataSet, State, DeuterationType, Peptides
import narwhals as nw
import pandas as pd

from pyhdx.process import correct_d_uptake, apply_control


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
    df = df.with_columns((nw.col("end") + 1).alias("stop"))
    return df


def load_pyhdx_peptides(peptides: list[Peptides]) -> dict[str, pd.DataFrame]:
    """Load peptides from hdxms_datasets into pyhdx Peptides model"""

    found_types = set(p.deuteration_type for p in peptides)
    if len(found_types) > len(peptides):
        raise ValueError("Dataset contains multiple peptide with the same deuteration type")

    peptide_types = {
        DeuterationType.partially_deuterated: "experiment",
        DeuterationType.fully_deuterated: "fd_control",
        DeuterationType.non_deuterated: "nd_control",
    }

    output = {}
    for dtype, attr_name in peptide_types.items():
        p = get_peptides_by_type(peptides, dtype)

        if p is None:
            continue

        peptide_df = adapt_for_pyhdx(p.load()).to_pandas()
        output[attr_name] = peptide_df

    return output


def parse_dataset(dataset: HDXDataSet, drop_first: int) -> list[tuple[pd.DataFrame, dict]]:
    """Parse an HDXDataSet into a list of tuples of (peptides, metadata) for pyhdx"""

    output = []
    loaded_peptides = [load_pyhdx_peptides(state.peptides) for state in dataset.states]

    if all("fd_control" in lp for lp in loaded_peptides):
        print("pass")

    def find_fd_control(loaded_peptides) -> tuple[int, pd.DataFrame]:
        for i, lp in enumerate(loaded_peptides):
            if "fd_control" in lp:
                return i, lp["fd_control"]
        raise ValueError("No FD control found")

    idx, fd_control = find_fd_control(loaded_peptides)

    import warnings

    warnings.warn(
        f"Not all states have FD control, using FD control from state index {idx}: {dataset.states[idx].name}"
    )

    # add the fd_control to all loaded peptides that lack it
    for lp in loaded_peptides:
        if "fd_control" not in lp:
            lp["fd_control"] = fd_control

    for state, peptides in zip(dataset.states, loaded_peptides):
        pd_peptides = get_peptides_by_type(state.peptides, DeuterationType.partially_deuterated)
        assert pd_peptides is not None  # this never happens due to previous checks

        peptides_merge = apply_control(**peptides)  # type: ignore
        peptides_corrected = correct_d_uptake(
            peptides_merge,
            drop_first=drop_first,
            d_percentage=pd_peptides.d_percentage or 100.0,
        )

        metadata = {
            **state_kwargs(state),
            **peptides_kwargs(pd_peptides),
        }

        output.append((peptides_corrected, metadata))

    return output
