!#/bin/bash

### Example start
curl https://ztf.snad.space/dr4/csv/633207400004730 | # Get some ZTF data
tail +2 | # chomp CSV header
awk -F, '{print $3"\t"$4"\t"$5}' | # print needed columns and change separator to tab
dmdt \
  --max-abs-dm=1.5 --height=64 \
  --min-lgdt=0 --max-lgdt=2 --width=96 \
  --smear --approx-smearing \
  --norm=lgdt --norm=max \
  --output=example.png
### Example end
