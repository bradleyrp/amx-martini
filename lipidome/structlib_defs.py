
defs_pts_lipids = [
	[0,.5,0,0,.5,0,0,.5,0,0,0,0,0,0,1,1,1,1,1,1],
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
	[10,9,9,8,8,7,6,6,5,4,3,2,1,0,5,4,3,2,1,0]]

#---! note that the "insane" definitions vs latest martini lipidome use CN0 vs CN0
#---! ...I was using only a portion of this code for less than five minutes before I detected this sloppy bug

defs_lipids = """
DTPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - 
DLPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - 
DPPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DBPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
POPC -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DOPC -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - 
DAPC -   -   -  NC3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - 
DIPC -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - 
DGPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - 
DNPC -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B
DTPE -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - 
DLPE -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - 
DPPE -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DBPE -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
POPE -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DOPE -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - 
POPG -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DOPG -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - 
POPS -   -   -  CNO  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - 
DOPS -   -   -  CNO  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - 
DPSM -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B  -   - 
DBSM -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B  - 
BNSM -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B C6B
OPPG -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B D2B C3B C4B  -   - 
JPPG -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - 
JFPG -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - 
GMO  -   -   -   -   -   -  GL1 GL2 C1A C2A D3A C4A C5A  -   -   -   -   -   -   - 
DHPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - 
DMPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - 
DSPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
POPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - 
DOPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - 
DUPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - 
DEPC.o -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B
DHPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - 
DLPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - 
DMPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - 
DSPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
POPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - 
DOPE.o -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - 
PPCS.o -   -   -  NC3  -  PO4 AM1 AM2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - 
DOPG.o -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - 
POPG.o -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - 
DOPS.o -   -   -  CNO  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - 
POPS.o -   -   -  CNO  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - 
CPG.o -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B  -   - 
PPG.o -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - 
PPT.o -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - 
DSMG.o -   -   -  C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
DSDG.oC61 C41 C11 C62 C42 C12 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
DSSQ.o -   -   S6 C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - 
"""

defs_pts_ptdins = [
	[.5,.5,0,0,1,.5,0,0,.5,0,0,0,0,0,0,1,1,1,1,1,1],
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
	[8,9,9,7,10,10,10,6,6,5,4,3,2,1,0,5,4,3,2,1,0]]

defs_ptdins = """
DPPI C1   C2   C3    PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
POPI C1   C2   C3    PO4   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
PIPI C1   C2   C3    PO4   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
PAPI C1   C2   C3    PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - 
PUPI C1   C2   C3    PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A   -   C1B  C2B  C3B  C4B   -    - 
POP1 C1   C2   C3    PO4  P1   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
POP2.o2 C1   C2   C3    PO4  P1  P2   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
POP2 C1   C2   C3    PO4  P1  P2   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
POP3.o2 C1   C2   C3    PO4  P1  P2  P3  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
POP3 C1   C2   C3    PO4  P1  P2  P3  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - 
PI.o C1   C2   C3    PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - 
PI34.o C1   C2   C3    PO4 PO1 PO2   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - 
"""

defs_ptdins_changes = """
changed ptdins lipid atoms named CP to PO4 per martini_v2.0_lipids_all_201506.itp
moved POP2 to POP2.o2 and got the topology manually. only change is C2A->D2A and D3A->C3A
ditto POP3
"""

defs_pts_sterols =[
	[0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,5.3,4.5,3.9,3.3,3,2.6,1.4,0,0,0,0,0]]

defs_sterols = """
CHOL  -   -   -   -   -   -  ROH R1  R2  R3  R4  R5  C1  C2   -   -   -   - 
ERGO - - - - - - ROH R1 R2 R3 R4 R5 C1 C2 - - - - 
"""