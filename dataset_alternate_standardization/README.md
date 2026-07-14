# Alternate standardization of FartDB

This folder contains a **robustness analysis** of the dataset curation and model evaluation reported in the FART paper [1].
The analysis was produced in response to a Matters Arising by Burns *et al.*.

The experiments here answer a single question:

> **If we re-curate FartDB with a much more aggressive standardization protocol — one that merges salt forms, protonation states, tautomers and stereoisomers into a single entry — do the dataset and the reported model performance change in any meaningful way?**

Short answer: **no.** The aggressive protocol removes ~1.4 % of the dataset and < 1 % of the test set, and every performance metric is essentially unchanged (see [Results](#results)).

---

## Why FartDB is standardized the way it is

Standardization for a *taste* dataset is not the same problem as standardization for a generic cheminformatics pipeline, because several transformations that are normally considered "clean-up" actually **erase information that determines taste**. Our guiding rule during curation was therefore:

> *Normalize the string representation of a molecule, but never merge two entries that a taste receptor would tell apart.*

**Applied in FartDB — pure representation clean-up (does not change taste identity):**

| Step | Rationale |
|---|---|
| Remove invalid SMILES | the string does not encode a real molecular graph |
| Remove multi-fragment entries | isolate the tastant from co-formulated solvents/salts |
| Remove net-charged entries without a counter-ion | a bare charge is an incomplete, ambiguous structure |
| Canonicalize SMILES (RDKit `StandardizeSmiles`) | one canonical string per molecule |
| Remove exact duplicates (same SMILES **and** taste) | the identical molecule–taste pair already exists |
| Remove MW > 2000 Da | exceeds the ChemBERTa context window |

The full core pipeline is in [`../dataset/scripts/FART_Data_Curation.ipynb`](../dataset/scripts/FART_Data_Curation.ipynb).

**Applied only by the extended protocol — each collapses taste-relevant chemistry:**

| Operation | What it would destroy | Concrete example |
|---|---|---|
| Strip stereochemistry | enantiomers/diastereomers can taste differently | L-glutamate is umami, D-glutamate is not [2]; the two enantiomers of aspartame are sweet vs. bitter [3] |
| Canonicalize tautomers | distinct tautomers can be distinct tastants | fructose: the β-pyranose ring tautomer is much sweeter than its furanose forms [4] |
| Salt stripping / reionization | protonation and salt form change perceived taste | glutamic acid: umami as its Na⁺ salt (MSG), sour when protonated [5] |
| InChIKey-based de-duplication | merges all of the above into a single record | removes the very distinctions the task depends on |

We apply these four operations only to *test robustness*, not to build the released dataset.

---

### A note on tautomers

Tautomer canonicalization is the subtlest of the excluded steps and the one most often assumed to be harmless. Taste can discriminate between forms of a molecule that differ only slightly, so different tautomers of one compound can taste different.

Fructose exists as an equilibrium of ring (pyranose/furanose) tautomers, and the β-D-fructopyranose form is markedly sweeter than the furanose forms [4]. As heating shifts the equilibrium toward the furanose forms, fructose (and honey) tastes less sweet warm than cold [4, 6]. In direct human tasting, the anomers of mannose — which interconvert by mutarotation, through the open-chain tautomer — differ categorically: α-D-mannose is sweet (sucrose-like) and β-D-mannose is bitter (quinine-like) [7].

For FartDB the risk is concrete. α-Amino acids undergo amido–imidol tautomerism, and RDKit's tautomer canonicalization removes sp3 stereochemistry involved in tautomerism *by default* [8]. As a result the step maps both L- and D-glutamic acid to the same achiral structure (`NC(CCC(=O)O)C(=O)O`), erasing the distinction between L-glutamate, which elicits umami, and D-glutamate, which does not [2].

Finally, a "canonical tautomer" is a bookkeeping device, not a physical claim: the RDKit documentation notes it is very likely *not* the most stable tautomer under any given conditions, only a reproducible one [8]. Sitzmann et al. [9] keep a canonical *identifier* while explicitly retaining the original tautomeric forms, and Baker et al. [10] show that the tautomer one standardizes to measurably affects downstream predictive models. Tautomers do interconvert, so *light* normalization of trivial drawing artifacts is reasonable; the point is only that *aggressive* canonicalization to a single heuristic tautomer is inappropriate for a taste dataset.

---

## The extended standardization procedure

`alternate_standardization.ipynb` re-processes the curated dataset
(`../dataset/fart_curated.csv`, 15,025 entries) with the following additional steps.

**1. Remove contradictory-label entries.**
Entries whose `Original Labels` field carries an explicit contradiction (labels containing `"not"` / `"contra"`) are flagged (17 candidates). After manual inspection, 5 unambiguously contradictory entries are removed.

**2. Extended per-molecule standardization.** For every molecule:

```python
from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

params = rdMolStandardize.CleanupParameters()
params.maxTransforms = 100000
params.maxTautomers  = 10000

salt_remover = SaltRemover.SaltRemover()
uncharger    = rdMolStandardize.Uncharger(canonicalOrder=True)

def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = salt_remover.StripMol(mol, dontRemoveEverything=True)   # salt stripping
    mol = uncharger.uncharge(mol)                                 # reionization / neutralization
    mol = rdMolStandardize.CanonicalTautomer(mol, params)         # tautomer canonicalization
    smiles_flat = Chem.MolToSmiles(mol, isomericSmiles=False)     # remove stereochemistry
    return rdMolStandardize.StandardizeSmiles(smiles_flat)        # final normalization
```

| Sub-step | Effect |
|---|---|
| `SaltRemover.StripMol` | removes salt/counter-ion fragments, keeping the parent tastant |
| `Uncharger.uncharge` | neutralizes formal charges (collapses protonation states / zwitterions) |
| `CanonicalTautomer` | maps every molecule to a single canonical tautomer |
| `MolToSmiles(isomericSmiles=False)` | discards stereochemistry (enantiomers/diastereomers merge) |
| `StandardizeSmiles` | final RDKit standardization |

**3. InChIKey-based de-duplication.**
Duplicates are now removed on `(InChIKey, taste)` instead of `(SMILES, taste)`.
The InChIKey is computed from the flattened, neutralized, tautomer-canonicalized structure, so it merges any residual salt protonation/tautomer/stereo variants that survived the string-level steps.
A `multi_taste` flag is recomputed by grouping on InChIKey.

**Output:** `fart_alternate_standardization.csv` — **14,810 entries**
(213 fewer than the 15,025 of FartDB, i.e. **−1.43 %**).

**4. Test-set cleaning.**
The held-out test set (`../dataset/splits/fart_test.csv`, 2,254 molecules) is put through the *same* pipeline and additionally filtered so that

- removed contradictory molecules are dropped,
- multi-taste molecules are removed, and
- only molecules that still exist in the extended training universe are kept.

This removes any near-duplicate that the more aggressive standardization would
expose as train/test overlap.

**Output:** `test_alternate_standardization.csv` — **2,237 molecules** (17 fewer than 2,254, i.e. **−0.75 %**).

---

## Files overview

| File | Description |
|---|---|
| `alternate_standardization.ipynb` | Builds the extended dataset and cleaned test set; produces the class-balance figure. |
| `results_with_alternate_test_set.ipynb` | Re-evaluates the trained FART / FART-augmented models on the cleaned test set (confusion matrix, per-class metrics, ROC/AUROC). |
| `fart_alternate_standardization.csv` | Extended-standardized dataset |
| `test_alternate_standardization.csv` | Extended-standardized, leakage-cleaned test set |


---

## References

1. Zimmermann, Y., Sieben, L., Seng, H., Pestlin, P. & Görlich, F. [A chemical language model for molecular taste prediction](https://doi.org/10.1038/s41538-025-00474-z). *npj Sci. Food* **9**, 122 (2025).
2. Nelson, G. et al. [An amino-acid taste receptor](https://doi.org/10.1038/nature726). *Nature* **416**, 199–202 (2002).
3. Mazur, R. H., Schlatter, J. M. & Goldkamp, A. H. [Structure–taste relationships of some dipeptides](https://doi.org/10.1021/ja01038a046). *J. Am. Chem. Soc.* **91**, 2684–2691 (1969).
4. Shallenberger, R. S. [Sugar structure and taste](https://doi.org/10.1021/ba-1971-0117.ch015). In *Carbohydrates in Solution* (Advances in Chemistry, Vol. 117), 256–263 (1973).
5. Kurihara, K. [Glutamate: from discovery as a food flavor to role as a basic taste (umami)](https://doi.org/10.3945/ajcn.2009.27462D). *Am. J. Clin. Nutr.* **90**, 719S–722S (2009).
6. Hyvönen, L., Varo, P. & Koivistoinen, P. Tautomeric equilibria of D-glucose and D-fructose: [polarimetric measurements](https://doi.org/10.1111/j.1365-2621.1977.tb12570.x) (*J. Food Sci.* **42**, 652–653, 1977) and [NMR spectrometric measurements](https://doi.org/10.1111/j.1365-2621.1977.tb12572.x) (*J. Food Sci.* **42**, 657–659, 1977).
7. Stewart, R. A., Carrico, C. K., Webster, R. L. & Steinhardt, R. G. Jr. [Physicochemical stereospecificity in taste perception of α-D-mannose and β-D-mannose](https://doi.org/10.1038/234220a0). *Nature* **234**, 220 (1971).
8. Landrum, G. et al. RDKit: Open-source cheminformatics. RDKit documentation — [`rdMolStandardize.TautomerEnumerator`](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html#rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator) (`Canonicalize`, `SetRemoveSp3Stereo`). https://www.rdkit.org (accessed July 2026).
9. Sitzmann, M., Ihlenfeldt, W.-D. & Nicklaus, M. C. [Tautomerism in large databases](https://doi.org/10.1007/s10822-010-9346-4). *J. Comput.-Aided Mol. Des.* **24**, 521–551 (2010).
10. Baker, C. M. et al. [Tautomer standardization in chemical databases: deriving business rules from quantum chemistry](https://doi.org/10.1021/acs.jcim.0c00232). *J. Chem. Inf. Model.* **60**, 3781–3791 (2020).