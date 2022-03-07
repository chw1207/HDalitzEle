import math

def deltaPhi(phi1, phi2):
    dPhi = phi1 - phi2
    if (dPhi > math.pi):
        dPhi -= 2. * math.pi
    if (dPhi < -math.pi):
        dPhi += 2. * math.pi
    return dPhi

def deltaR(eta1, phi1, eta2, phi2):
    dEta = eta1 - eta2
    dPhi = deltaPhi(phi1, phi2)
    return math.sqrt(dEta * dEta + dPhi * dPhi)

def deltaEta(eta1, eta2):
    return eta1 - eta2

def twoBodyPt(pt1, phi1, pt2, phi2):
    px1 = pt1 * math.cos(phi1)
    py1 = pt1 * math.sin(phi1)
    px2 = pt2 * math.cos(phi2)
    py2 = pt2 * math.sin(phi2)
    return math.sqrt((px1+px2)**2 + (py1+py2)**2)