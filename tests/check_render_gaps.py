#!/usr/bin/env python3
"""Check two clear ring-gap pixels in the six-track synthetic smoke fixture.

At y=1500, x=707 and 808 lie midway in the 0.02r gaps next to the two
histograms on scaffold_T1, whose TE values are all zero. Both native and
SVG-derived PNG should be white there. This detects an observed rendering
mismatch; it does not establish general visual fidelity.
"""
from pathlib import Path
import json
import struct
import sys
import zlib

POINTS = [(707, 1500), (808, 1500)]


def pixels(path):
    data=Path(path).read_bytes()
    assert data[:8]==b'\x89PNG\r\n\x1a\n'
    w,h,depth,color,compression,filtering,interlace=struct.unpack('>IIBBBBB',data[16:29])
    assert (w,h)==(3000,3000) and depth==8 and color in (2,6) and interlace==0
    packed=bytearray(); offset=8
    while offset<len(data):
        size=struct.unpack('>I',data[offset:offset+4])[0]; kind=data[offset+4:offset+8]
        if kind==b'IDAT': packed.extend(data[offset+8:offset+8+size])
        offset+=size+12
    raw=zlib.decompress(packed); channels=3 if color==2 else 4; stride=w*channels
    previous=bytearray(stride); offset=0; result={}
    def paeth(a,b,c):
        p=a+b-c; errors=(abs(p-a),abs(p-b),abs(p-c))
        return (a,b,c)[errors.index(min(errors))]
    for y in range(max(y for _,y in POINTS)+1):
        kind=raw[offset]; row=bytearray(raw[offset+1:offset+1+stride]); offset+=stride+1
        assert 0<=kind<=4
        for i in range(stride):
            a=row[i-channels] if i>=channels else 0; b=previous[i]; c=previous[i-channels] if i>=channels else 0
            prediction=(0,a,b,(a+b)//2,paeth(a,b,c))[kind]
            row[i]=(row[i]+prediction)&255
        for x,py in POINTS:
            if py==y: result[f'{x},{y}']=list(row[x*channels:x*channels+3])
        previous=row
    return result


if __name__=='__main__':
    if len(sys.argv)!=3: raise SystemExit('Usage: check_render_gaps.py NATIVE.png SVG_DERIVED.png')
    results={name:pixels(path) for name,path in zip(('native','svg'),sys.argv[1:])}
    print(json.dumps(results))
    if not all(min(rgb)>=250 for samples in results.values() for rgb in samples.values()):
        raise SystemExit('FAIL: expected white synthetic ring gaps; native/SVG rendering fidelity needs review')
