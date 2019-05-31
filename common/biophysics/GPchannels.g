// Load ion channel description files
float xmin = -0.2
float xmax = 0.1
int xdivs = 300
float dx = (xmax - xmin)/xdivs
float x
int i

include ../common/biophysics/NaF.g
include ../common/biophysics/NaP.g
include ../common/biophysics/Kv2.g
include ../common/biophysics/Kv3.g
include ../common/biophysics/Kv4.g
include ../common/biophysics/KCNQ.g
include ../common/biophysics/SK.g
include ../common/biophysics/CaHVA.g
include ../common/biophysics/CaBuff.g
include ../common/biophysics/HCN.g
