from __future__ import annotations

import time
import timeit
from typing import Dict, List, Optional, Sequence, Set, Tuple, cast

from rdkit import Chem
from rdkit.Chem import AllChem

from rdcanon.main import canon_reaction_smarts, canon_smarts, random_smarts


def compare_reaction_outputs(
    reactant_objs_in: List[Chem.Mol],
    template_list: List[Chem.rdChemReactions.ChemicalReaction],
    canon_template_list: List[Chem.rdChemReactions.ChemicalReaction],
) -> Tuple[int, int]:
    correct, incorrect, failed = 0, 0, 0
    ordered_arrs: List[Tuple[Chem.Mol, ...]] = []
    for i2, k in enumerate(template_list):
        ordered_noncanon = (reactant_objs_in[i2],)
        ordered_canon = (reactant_objs_in[i2],)

        p = k.RunReactants(ordered_noncanon)
        p2 = canon_template_list[i2].RunReactants(ordered_canon)

        p1s = []
        p1sms = []
        for l in p:
            for ll in l:
                p1s.append(ll)
                sm_out = Chem.MolToSmiles(ll, isomericSmiles=False)
                p1sms.append(sm_out)

        p2s = []
        p2sms = []
        for l in p2:
            for ll in l:
                p2s.append(ll)
                sm_out = Chem.MolToSmiles(ll, isomericSmiles=False)
                p2sms.append(sm_out)

        all_hit = True
        for ppp1, ppp2 in zip(sorted(p1sms), sorted(p2sms)):
            if ppp1 != ppp2:
                all_hit = False
                break

        if all_hit:
            correct = correct + 1
        else:
            pass
            incorrect = incorrect + 1
    return correct, incorrect


def compare_products(reaction_template: str, reactants_in: str) -> bool:
    canon_rxn1 = canon_reaction_smarts(reaction_template, True, "drugbank", True)
    rxn = AllChem.ReactionFromSmarts(reaction_template)  # type: ignore[attr-defined]
    rxn_canon = AllChem.ReactionFromSmarts(canon_rxn1)  # type: ignore[attr-defined]

    reactants = Chem.MolFromSmiles(reactants_in)
    p = rxn.RunReactants((reactants,))
    p2 = rxn_canon.RunReactants((reactants,))

    if len(p) != len(p2):
        return False

    p1s = []
    p1sms = []
    for l in p:
        for ll in l:
            Chem.SanitizeMol(ll)
            p1s.append(ll)
            sm_out = Chem.MolToSmiles(ll, isomericSmiles=False)
            sm_canon = Chem.CanonSmiles(sm_out)
            p1sms.append(sm_canon)

    p2s = []
    p2sms = []
    for l in p2:
        for ll in l:
            Chem.SanitizeMol(ll)
            p2s.append(ll)
            sm_out = Chem.MolToSmiles(ll, isomericSmiles=False)
            sm_canon = Chem.CanonSmiles(sm_out)
            p2sms.append(sm_canon)

    all_hit = True
    for ppp1, ppp2 in zip(sorted(p1sms), sorted(p2sms)):
        if ppp1 != ppp2:
            all_hit = False
            break

    if all_hit:
        return True
    else:
        return False


def find_n_matches(
    smarts_library: List[Chem.Mol],
    target_library: Sequence[Optional[Chem.Mol]],
    n: int,
) -> Dict[str, Dict[str, List[str]]]:
    out_data: Dict[str, Dict[str, List[str]]] = {}

    for smol in smarts_library:
        sm = Chem.MolToSmarts(smol)
        if sm not in out_data:
            out_data[sm] = {
                "query_smarts": [],
                "matching_substrate_smiles": [],
                "non_matching_substrate_smiles": [],
                "random_substrate_smiles": [],
            }

        match_hit = False
        non_match_hit = False
        match_smiles = ""
        non_match_smiles = ""
        tot = 0
        for r in target_library:
            if r is None:
                continue

            if r.HasSubstructMatch(smol):
                match_hit = True
                match_smiles = Chem.MolToSmiles(r)
            else:
                non_match_hit = True
                non_match_smiles = Chem.MolToSmiles(r)

            if match_hit and non_match_hit:
                out_data[sm]["query_smarts"].append(sm)
                out_data[sm]["matching_substrate_smiles"].append(match_smiles)
                out_data[sm]["non_matching_substrate_smiles"].append(non_match_smiles)
                match_hit = False
                non_match_hit = False
                match_smiles = ""
                non_match_smiles = ""
                tot = tot + 1
            if tot == n:
                break

    return out_data


def run_against_library(
    smarts_library: List[str],
    target_library: List[str],
    n: int,
    emb: str = "drugbank",
) -> Tuple[Dict[str, Dict[str, List[str]]], Dict[str, Dict[str, List[str]]]]:
    noncanon_template_obj = [
        Chem.MolFromSmarts(template_smarts) for template_smarts in smarts_library
    ]
    mols = [Chem.MolFromSmiles(substrate_smiles) for substrate_smiles in target_library]
    noncanon_output = find_n_matches(noncanon_template_obj, mols, n)

    canon_template_obj = [
        Chem.MolFromSmarts(canon_smarts(template_smarts, mapping=True, embedding=emb))
        for template_smarts in smarts_library
    ]

    canon_output = find_n_matches(canon_template_obj, mols, n)

    return noncanon_output, canon_output


def compare_product_sets(
    noncanon_output: Dict[str, Dict[str, List[str]]],
    canon_output: Dict[str, Dict[str, List[str]]],
) -> bool:
    for noncanon_hit, canon_hit in zip(
        noncanon_output,
        canon_output,
    ):
        if (
            noncanon_output[noncanon_hit]["matching_substrate_smiles"]
            != canon_output[canon_hit]["matching_substrate_smiles"]
        ):
            print(noncanon_hit, canon_hit)
            return False

        if (
            noncanon_output[noncanon_hit]["non_matching_substrate_smiles"]
            != canon_output[canon_hit]["non_matching_substrate_smiles"]
        ):
            print(noncanon_hit, canon_hit)
            return False

    return True


def run_random_permutations(in_smarts: str, n_perms: int = 100) -> bool:
    all_random: Set[str] = set()
    for i in range(n_perms):
        sm = random_smarts(in_smarts)
        # print(sm)
        all_random.add(sm)

    canon_out: List[str] = []
    for r in all_random:
        canon_out.append(cast(str, canon_smarts(r)))

    return len(set(canon_out)) == 1


def compare_substrate_datasets(
    template_smarts_dataset: List[Chem.Mol],
    substrate_smiles_dataset: List[Chem.Mol],
) -> None:
    for template_smarts, substrate_smiles in zip(
        template_smarts_dataset, substrate_smiles_dataset
    ):
        substrate_smiles.HasSubstructMatch(template_smarts)


def time_compare_substruct_match(
    template_smarts_dataset: List[str],
    substrate_smiles_dataset: List[str],
    embeddings: List[str] = ["drugbank"],
    iters: int = 10,
    v: bool = False,
) -> List[float]:
    noncanon_template_obj = [
        Chem.MolFromSmarts(template_smarts)
        for template_smarts in template_smarts_dataset
    ]
    mols = [
        Chem.MolFromSmiles(substrate_smiles)
        for substrate_smiles in substrate_smiles_dataset
    ]

    t1 = timeit.timeit(
        lambda: compare_substrate_datasets(noncanon_template_obj, mols),
        number=iters,
        timer=time.process_time,
    )

    ts = [t1]
    embeds = ["rdchiral"]
    for idx, emb in enumerate(embeddings):
        canon_template_obj = [
            Chem.MolFromSmarts(
                canon_smarts(template_smarts, mapping=True, embedding=emb)
            )
            for template_smarts in template_smarts_dataset
        ]

        t2 = timeit.timeit(
            lambda: compare_substrate_datasets(canon_template_obj, mols),
            number=iters,
            timer=time.process_time,
        )

        ts.append(t2)
        embeds.append("embed_" + str(idx + 2))

    if v:
        for t, e in zip(ts, embeds):
            print(e, t)
    return ts
