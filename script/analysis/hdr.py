import hdf5_to_dict as io
import sys

if len(sys.argv) != 2:
  print('ERROR format is')
  print('  python hdr.py [filename]')
  sys.exit()

fnam = sys.argv[1]

hdr = io.load_hdr(fnam)

#print(hdr.keys())

def print_val(vnam):
  if vnam in hdr.keys():
    print(vnam,'=','%e' % hdr[vnam])

print_val('N1')
print_val('N2')
print_val('N3')
print_val('a')
print_val('M_unit')
print_val('B_unit')
print_val('L_unit')
print_val('Ne_unit')

