
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
gcd = graph_codec_rdkit()
x = "[!#8!#15!#16!#17!#35!+!X4:1]~[*:2]"
print(x)
sma = gcd.smarts_decode(x)
print(gcd.smarts_encode(sma))

