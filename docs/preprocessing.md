# Preprocessing

Raw peak lists are preserved unchanged in the normalized dataset. Parser normalization is limited
to typed numeric peak representation and metadata-key normalization; it does not rewrite source
intensities or peak locations.

## Pilot baseline preprocessing

The fixed pilot EI-MS baselines use a deliberately simple, frozen representation applied within
the modeling layer only:

1. round/assign peaks to nominal integer m/z bins from 0 through 700;
2. sum intensity if multiple source peaks land in the same nominal bin;
3. apply a square-root intensity transform;
4. L2-normalize each spectrum vector.

No peak filtering, learned embedding, alignment optimization, or fold-specific preprocessing
hyperparameter search is performed. The spectral 1-NN baseline uses cosine similarity on these
vectors; logistic regression, linear SVM, and random forest use the same fixed representation.

These choices define the pilot baselines only. Future benchmark models may use other documented
representations, but they must not alter the preserved raw peak list and must fit any learned
preprocessing using training folds only.
