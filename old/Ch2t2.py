import string

original = "g fmnc wms bgblr rpylqjyrc gr zw fylb. rfyrq ufyr amknsrcpq ypc " \
    "dmp. bmgle gr gl zw fylb gq glcddgagclr ylb rfyr'q ufw rfgq " \
    "rcvr gq qm jmle. sqgle qrpgle.kyicrpylq() gq pcamkkclbcb. lmu " \
    "ynnjw ml rfc spj.http://www.pythonchallenge.com/pc/def/map.html"\
    "g fmnc wms bgblr rpylqjyrc gr zw fylb. rfyrq ufyr amknsrcpq ypc dmp. bmgle gr gl zw fylb gq glcddgagclr ylb rfyr'q ufw rfgq rcvr gq qm jmle. sqgle qrpgle.kyicrpylq() gq pcamkkclbcb. lmu ynnjw ml rfc dgjclykc gl rfc spj."\
    "http://www.pythonchallenge.com/pc/def/map.html"

table = string.maketrans(
    "abcdefghijklmnopqrstuvwxyz", "cdefghijklmnopqrstuvwxyzab"
)

print original.translate(table) 