import pandas as pd
from clockcheck import LinearClock
from biolearn.data_library import DataLibrary

print("Downloading GSE64495 (large file, first run is slow, it caches)...")
data = DataLibrary().get("GSE64495").load()
print("Loaded. Metadata columns:", list(data.metadata.columns))

mine = LinearClock.load("Horvath1").predict(data.dnam, on_missing="drop")
theirs = pd.to_numeric(data.metadata["dna_methylation_age"], errors="coerce")

comp = mine.to_frame("clockcheck").join(theirs.rename("published")).dropna()
print("n samples compared:", len(comp))
print("correlation:", round(comp.corr().iloc[0, 1], 4))
print("median abs diff (years):", round((comp.clockcheck - comp.published).abs().median(), 2))
print(comp.head(10))
