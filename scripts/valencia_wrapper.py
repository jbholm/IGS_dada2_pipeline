#!/usr/local/packages/python-3.5.2/bin/python
import sys, os, subprocess
import pandas as pd

input = sys.argv[2]
temp = input + ".temp"
df = pd.read_csv(input, header=0, index_col=0, sep=',')
# Calculate the total read count for each sample and insert as 2nd column
df.insert(0, 'read_count', df.sum(axis=1))
df.to_csv(temp)

process = subprocess.run([os.path.join(sys.path[0], "Valencia_v1.py"), sys.argv[1], temp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# if Valencia exited 1, print Valencia's stderr and raise CalledProcessError to parent
if process.returncode:
    print(process.stderr.decode(sys.stdout.encoding))
    process.check_returncode()  #Raises CalledProcessError which itself contains stderr
else:
    os.remove(temp)
    sys.exit(process.returncode)
