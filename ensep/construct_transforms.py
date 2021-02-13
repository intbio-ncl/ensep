"""
This script takes the full EnSeP MSGPack files and constructs smaller
files for individual transforms.
"""

from itertools import chain
from multiprocessing import Pool, cpu_count

import molvs
import msgpack
import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

standardizer = molvs.standardize.Standardizer(tautomer_transforms=[])

with open("ensep.msgpack", "rb") as f:
    d = msgpack.unpackb(f.read(), raw=False, use_list=False)


def standardize(mol):
    try:
        sp_mol = standardizer.super_parent(mol)
        sp_smiles = Chem.MolToSmiles(sp_mol)
        return sp_smiles
    except ValueError:
        return None


def test(a_mol, b_smiles, transform):
    a_prods = transform.RunReactants([a_mol])
    if not a_prods:
        return False

    a_prods = [standardize(i) for i in chain(*a_prods)]

    if not any(["[*]" in i for i in a_prods]):
        if b_smiles in a_prods:
            return True
        return False

    else:
        qp = Chem.AdjustQueryParameters()
        qp.makeDummiesQueries = True
        qp.adjustDegree = True
        qp.adjustDegreeFlags = Chem.ADJUST_IGNOREDUMMIES

        a_prods = [Chem.AddHs(i) for i in a_prods]
        a_prods = [Chem.AdjustQueryProperties(i, qp) for i in a_prods]

        b_mol = Chem.MolFromSmiles(b_smiles)
        b_mol = Chem.AddHs(b_mol)

        if any([b_mol.HasSubstructMatch(i) for i in a_prods]):
            return True

        return False


def run_single(tup):
    a, b, tfm = tup

    a_mol = Chem.MolFromSmiles(a)
    a_mol = Chem.AddHs(a_mol)

    b_mol = Chem.MolFromSmiles(b)
    b_mol = Chem.AddHs(b_mol)

    if test(a_mol, b, tfm):
        return (a, b)
    elif test(b_mol, a, tfm):
        return (b, a)
    
    return None

def run(transform):
    global d
    tfm = AllChem.ReactionFromSmarts(transform)

    results = {}

    total = len(d)

    jobs = ((a, b, tfm) for a, b in d.keys())

    P = Pool(cpu_count())
    for i in tqdm.tqdm(P.imap_unordered(run_single, jobs), total=total):
        if i:
            results[i] = d[tuple(sorted(i))]

    return results


tfms = {
    "ecrc": "[c:1]1([O:7])[c:2]([O:8])[c:3][c:4][c:5][c:6]1>>[C:1](-[*:7])(=[O])[C:2]([*:8])=[C:3]-[C:4]=[C:5]-[C:6](=[O])",
    "icrc": "[c:1]1([O:7][H])[c:2][c:3][c:4][c:5][c:6]([O:8][H])1>>[C:1](=[O:7])(-[O])[C:2]=[C:3]-[C:4]=[C:5]-[C:6](=[O:8])(-[O])",
    "phenol_hydrox": "[c:1][c:2][c:3]([O:7])[c:4]([H:8])[c:5][c:6]>>[c:1][c:2][c:3]([O:7])[c:4]([O:8])[c:5][c:6]",
    "o2a": "[C:2]([C:1])([C:3])([H:4])-O>>[C:2]([*:1])([*:3])=O",
    "oxidation_CC": "[C:2]([H:1])-[C:3]([H:4])>>[C:2]=[C:3]",
    "ee2_cleave_1": "[*:1]#[*:2]>>[*:1](=[#8])[#8]",
    "ee2_cleave_2": "[*:1]#[*:2]>>[*:1](-[*:2]([H])[H])=[O]",
}

for k, v in tfms.items():
    results = run(v)
    print(k, len(results))
    with open(f"2020_05_11_{k}.msgpack", "wb") as f:
        f.write(msgpack.packb(results))
