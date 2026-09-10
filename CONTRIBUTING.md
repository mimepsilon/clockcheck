# Contributing

The bar for this project is honesty, not volume. A good contribution makes a
clock's failure modes *more* visible, never less.

Ground rules:
- **Never commit fabricated coefficients or data.** Clock weights come from
  published sources; datasets are real (GEO/your own) or the labeled synthetic
  fixture. Synthetic output must stay stamped `NON-EVIDENTIAL`.
- **A new audit ships with a test** on the synthetic fixture that shows it *fires*
  when the failure is present and stays quiet when it isn't.
- **Caveats first.** Report layers lead with what could make a number wrong.

Dev setup:
```bash
pip install -e ".[dev,real]"
pytest -q
```
