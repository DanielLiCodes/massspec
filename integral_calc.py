
from gnpsdata.workflow_classicnetworking import get_clustersummary_dataframe
df = get_clustersummary_dataframe("b53687f90fcb4e4f8b4dbc3d5ab2a376")
import pymzml
import bisect

mzml_file = 'Hui_N1.mzML'
run = pymzml.run.Reader(mzml_file)

time_mz_i = []
for n, spec in enumerate(run):
  if spec.ms_level != 1:
    continue
  
  mz_i = spec.peaks('centroided')

  for i in range(1, len(mz_i)):
    mz_i[i][1] += mz_i[i - 1][1]

  mz = [-10000]
  i = [0]
  for mass, intensity in mz_i:
    mz.append(mass)
    i.append(intensity)

  time_mz_i.append((spec.scan_time_in_minutes(), (mz, i)))

result = {}
time_tolerance = 0.5 # minutes
mz_tolerance = 0.5 # Da
for index, row in df.iterrows():
  area = 0
  previous_time = 0
  previous_intensity = 0
  for time, (mz, i) in time_mz_i:
    if abs(time - row['RTMean'] / 60) > time_tolerance:
      previous_time = time
      previous_intensity = 0
      continue
    
    lower = bisect.bisect_left(mz, row['precursor mass'] - mz_tolerance) - 1
    upper = bisect.bisect_right(mz, row['precursor mass'] + mz_tolerance) - 1

    current_intensity = i[upper] - i[lower]
    
    area += (time - previous_time) * (current_intensity + previous_intensity) / 2
    previous_time = time
    previous_intensity = current_intensity

  result[row['cluster index']] = (row['precursor mass'], row['RTMean']/60,area)

print(len(result))

for key in result.keys():
  print(f'Cluster: {key:{4}}, Mass: {result[key][0] :8.3f}, Average time: {result[key][1] :.2f}, Area: {result[key][2]: 17.6f}')