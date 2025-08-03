import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

fileName = sys.argv[1]
title = sys.argv[2]
#functionName = sys.argv[2]
#numTests = sys.argv[3]
#equation = sys.argv[4]

fig, ax = plt.subplots(figsize=(16, 12))

df = pd.read_csv(fileName).drop(['Index'], axis=1)

styles = ["solid"]#, "dotted", "dashed", "dashdot"]
styleIndex = 0

for col in df.columns.tolist():
	ax.plot(df[col], label=col, linestyle=styles[styleIndex])
	styleIndex = (styleIndex + 1) % len(styles)

ax.legend()

plt.title(title)
#plt.title(functionName + ': ' + numTests + ' tests' + "\n" + equation)
plt.ylabel('Avg Abs Error')
plt.xlabel('Samples')

fig.axes[0].set_xlim(1)#, len(df))

fig.axes[0].set_xscale('log', base=10)
fig.axes[0].set_yscale('log', base=10)

fig.tight_layout()
fig.savefig(os.path.splitext(fileName)[0] + ".png", bbox_inches='tight')
