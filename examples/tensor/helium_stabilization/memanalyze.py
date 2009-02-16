#/usr/bin/env python

import sys

allocs = {}
allocate = 0
deallocate = 0

def inc(val):
	if val not in allocs:
		allocs[val] = 0
	allocs[val] += 1
	
def dec(val):
	if val not in allocs:
		allocs[val] = 0
	allocs[val] -= 1
	

f = file(sys.argv[1])
for l in f.readlines():
	if l.startswith("BLITZALLOCATE"):
		exec(l)
		allocate += BLITZALLOCATE
		inc(BLITZALLOCATE)
	elif l.startswith("BLITZDEALLOCATE"):
		exec(l)
		deallocate += BLITZDEALLOCATE
		dec(BLITZDEALLOCATE)
		
	
print "Allocated %i bytes" % (allocate)
print "Deallocated %i bytes" % (deallocate)
print "Leaked %i bytes" % (allocate - deallocate)

for amount, count in allocs.iteritems():
	if count != 0:
		print "    %i bytes alloc-dealloc %i times" % (amount, count)
