# plotting.py

#==========================================================

import numpy as np
import matplotlib.pyplot as plt

#==========================================================

# oneField
def oneField(f,title):

	plt.contourf(f); plt.colorbar()
	plt.title(title)
	plt.tight_layout()
	plt.show()

	return 

# twoFields
def twoFields(fields,titles):

	plt.figure(figsize=(15,6))	

	plt.subplot(121)
	plt.contourf(fields[0]); plt.colorbar()
	plt.title(titles[0])

	plt.subplot(122)
	plt.contourf(fields[1]); plt.colorbar()
	plt.title(titles[1])

	plt.tight_layout()
	plt.show()
