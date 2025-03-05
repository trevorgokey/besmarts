#!/usr/bin/env python3

import sys
import numpy as np

gpermol2au = 1836
au2amu = 1/1836.
amu2au = 1822.888184038791
kcal2hartree = 1/627.5096080306
ang2bohr = 1/0.529177249
bohr2angstrom = 0.529177249
hartree2kcalmol = 627.5096080306
# eig2wavenumber = 1.0 /(2*np.pi * 137 * 5.29e-9) * np.sqrt(au2amu)
eig2wavenumber = 33.3564095 / (2*np.pi) 
eig2thz = eig2wavenumber * 2.99792458e10 / 1e12 
# conversion = eig2wavenumber
au2cm = 5.29177249e-9
c = 137
conversion = 1.0 /(2*np.pi * c * au2cm)

def converteig( eigs):
    return np.sign( eigs) * np.sqrt( np.abs( eigs))*conversion

# swap_internals = not connect
debug = False

def printortho( mat):
    if not debug:
       return 
    for i in mat:
        for x in i:
            x = abs(x)
            if x > 1e-8:
                print("{:3.1f}".format(x), end="")
            else:
                print("Y", end="")
        print()

def printmat( mat, dec=4):
    if not debug:
       return 
    for row in mat:
        for val in row:
            print(("{:"+str(dec+3)+"."+str(dec)+"f} ").format( val), end="")
        print()


mass_table = {
    'C': 12.011,
    'H': 1.008,
    'N': 14.007,
    'O': 15.999,
    'P': 30.973,
    'S': 32.06,
    'F': 18.998,
    "Cl": 35.45,
    "B": 10.81,
    "I": 126.9,
    "Br": 79.04,
    "Na": 22.909,
    "Mg": 24.305,
    "Li": 6.9410,
    "K": 39.098,
    "Ca": 40.078,
    "He": 4.0026,
    "Ne": 20.180,
    "Ar": 39.948,
    "Si": 28.086,
    "Zn": 65.380,
    "Cd": 112.41,
    "Ni": 58.693,
    "Kr": 83.798,
}
for k, v in list(mass_table.items()):
    mass_table[k.upper()] = v

# TODO clean and move these GS to a better home
# The useful code uses QR so these are no longer needed here
def gram_schmidt(vectors):
    basis = []
    for v in vectors:
        w = v - np.sum( np.dot(v,b)*b  for b in basis )
        if (w > 1e-10).any():
            basis.append(w/np.linalg.norm(w))
    return np.array(basis)


def gs2(X, row_vecs=True, norm=True):
    if not row_vecs:
        X = X.T
    Y = X[0:6][:].copy()
    if debug:
        print("gs2 Y.shape")
        print(Y.shape)
    for i in range(6, X.shape[1]):
        v = np.zeros(X.shape[1])
        v[i] = 1.0
        proj = np.array([np.dot(v,Y[b])*Y[b] for b in range(i)]).sum(axis=0)
            #proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, v - proj))
        if norm:
            Y[-1] /= np.linalg.norm(Y[-1]) 
        #Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T


def gs(X, row_vecs=True, norm = True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T


def save_xyz(syms, xyz, outfd=sys.stdout, comment=""):
    outfd.write(str(len(syms)) + "\n")
    outfd.write(comment + "\n")
    fmt = "{:2s}     {:16.13f} {:16.13f} {:16.13f}\n"
    [outfd.write(fmt.format(sym, *pos)) for (sym, pos) in zip(syms, xyz)]


def mode_to_xyz(syms, node, mode, outfd=sys.stdout, comment="", norm=True, magnitude=10.0):

    # from Lee-Ping Wang's version
    if(norm):
        mode = mode / np.linalg.norm(mode)/np.sqrt(mode.shape[0])
    spac = np.linspace(0, 1, 21)

    # one complete period (up, down, down, up)
    disp = magnitude*np.concatenate((spac, spac[::-1][1:], -1*spac[1:], -1*spac[::-1][1:-1]))

    m = np.array([mass_table[sym] for sym in syms])[:,np.newaxis]
    mode = mode * bohr2angstrom#  * np.sqrt(m)
    weight = np.array([mass_table[sym] for sym in syms])
    #weight = np.ones(len(syms))
    weight /= weight.sum()
    nodeCOM =  (node * weight[:,np.newaxis]).sum(axis=0)
    for dx in disp:
        save_xyz(syms, node - nodeCOM + mode*dx, outfd=outfd, comment=comment)


def mode_COM_to_xyz(syms, node, mode, outfd=sys.stdout, comment="", norm=True, magnitude=0.1):

    # from Lee-Ping Wangs version
    if(norm):
        mode = mode / np.linalg.norm(mode)/np.sqrt(mode.shape[0])
    spac = np.linspace(0, 1, 21)

    # one complete period (up, down, down, up)
    disp = magnitude*np.concatenate((spac, spac[::-1][1:], -1*spac[1:], -1*spac[::-1][1:-1]))
    

    weight = np.array([mass_table[sym] for sym in syms])
    #weight = np.ones(len(syms))
    weight /= weight.sum()
    nodeCOM =  np.atleast_2d((node + node.mean(axis=0) * weight[:,np.newaxis]).sum(axis=0))
    modeCOM =  np.atleast_2d((mode * (weight[:,np.newaxis])).sum(axis=0))

    syms = ['X']

    for dx in disp:
        save_xyz(syms, nodeCOM + modeCOM*dx, outfd=outfd, comment=comment)


def angular_momentum( xyz, mass, mode):

    COM = center_of_mass( xyz, mass)
    xyz = xyz - COM

    L = np.cross( mass.reshape( -1, 1) * mode.reshape( -1, 3), xyz, axis=1)

    return np.linalg.norm(L.sum( axis=0))


def center_of_mass(xyz, mass):
    mass = mass / mass.sum()
    COM = ( xyz * mass[ :, np.newaxis]).sum( axis=0)
    return COM


def inertia_tensor(xyz, mass):
    #return inertia_tensor2(xyz,mass)
    I = np.zeros((3,3))
    COM = center_of_mass(xyz, mass)
    xyz = xyz - COM
    for (x, y, z), m in zip(xyz, mass):
        I[0,0] += m*(y*y + z*z)
        I[0,1] -= m*(x*y)
        I[0,2] -= m*(x*z)
        I[1,1] += m*(x*x + z*z)
        I[1,2] -= m*(y*z)
        I[2,2] += m*(x*x + y*y)

    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
    if debug:
        print("I:")
        print(I)
    return I


def trans_rot_modes(xyz, mass, X, rot=True):
    #X = np.eye(3)
    if rot:
        D = np.zeros((6, np.product(xyz.shape)))
    else:
        D = np.zeros((3, np.product(xyz.shape)))
    na, nd = xyz.shape
    if debug:
        print(D.shape, xyz.shape)
    sqrt_mass = np.sqrt(mass)
    xyz = xyz - center_of_mass(xyz, mass)
    if debug:
        print(np.tile([1, 0, 0], na).shape)
        print("sqrt_mass")
        print(sqrt_mass)
    D[0][:] = np.tile([1, 0, 0], na) * np.repeat(sqrt_mass, nd)
    D[1][:] = np.tile([0, 1, 0], na) * np.repeat(sqrt_mass, nd)
    D[2][:] = np.tile([0, 0, 1], na) * np.repeat(sqrt_mass, nd)
    if debug:
        print("D[2]")
        print(D[2])
    #P = np.array([np.dot(xyz, p_ax[i]) for i in range(nd)]).reshape(nd,na)
    P = np.dot(xyz, X)
    x,y,z = 0,1,2
    if debug:
        print("P.shape")
        print(P.shape)

    if not rot:
        return D

    D[3][:] = np.array([(P[i,y]*X[j,z] - P[i,z]*X[j,y]) * sqrt_mass[i] \
            for i in range(na) for j in range(nd)])
    D[4][:] = np.array([(P[i,z]*X[j,x] - P[i,x]*X[j,z]) * sqrt_mass[i] \
            for i in range(na) for j in range(nd)])
    D[5][:] = np.array([(P[i,x]*X[j,y] - P[i,y]*X[j,x]) * sqrt_mass[i] \
            for i in range(na) for j in range(nd)])

    D = np.array([d for d in D if np.linalg.norm(d) > 0.0])
    if debug:
        print("TR_MODES")
        for row in np.dot(D,D.T):
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()

        print("D[3,:15", D[3][:15])
        print("D SHAPE", D.shape)

    #dot = np.dot(xyz, X)
    #print("dot.shape")
    #print(dot.shape)
    #for i in range(3):
    #    D[i+3,:] = np.concatenate([np.cross(dot, X[j])[:,i] for j in range(3)])

    #D[3:6] = np.cross(np.dot(P,X),np.dot(P,X)).reshape(3,-1)
    #D[3] /= np.linalg.norm(D[3])
    #D[4] /= np.linalg.norm(D[3])
    #D[5] /= np.linalg.norm(D[3])
    return D

def hessian_transform_mass_weighted(hess, mass):

    # assumes hess is kcal/A/A using mass in amu

    mass = np.array(mass) * amu2au
    hess = np.array(hess)
    mass_mat = (np.dot(np.sqrt(mass).reshape(-1,1), np.sqrt(mass).reshape(1,-1)))

    hess_convert = 1.0 / mass_mat * kcal2hartree / ang2bohr**2

    hess *= hess_convert

    return hess


def hessian_modes( hess, syms, xyz, mass, mol_number, remove=0, stdtr=True, debug=False, verbose=False, return_DL=False):

    # xyz = xyz * ang2bohr
    xyz = np.array(xyz)
    convert = 1
    mass = np.array(mass)
    mass_per_atom = mass[:,0]
    hess = hessian_transform_mass_weighted(hess, mass)

    E = np.linalg.eigh(hess * convert)# 
    if debug:
        print("freqs before removing trans/rot")
        print(np.sign(E[0])*np.sqrt(np.abs(E[0]))*conversion)
    esrt = E[0].argsort()

    mol_name = "mol_" + str(mol_number)
    if verbose:
        for i, mode in enumerate(E[1]):
            with open(mol_name+".mode_"+str(i)+".xyz", 'w') as fid:
                mode_to_xyz(syms, xyz, mode.reshape(-1, 3), outfd=fid, 
                        comment=mol_name + " mode " + str(i))
        for i, mode in enumerate(E[1]):
            with open(mol_name+".COM.mode_"+str(i)+".xyz", 'w') as fid:
                mode_COM_to_xyz(syms, xyz, mode.reshape(-1, 3), outfd=fid, 
                        comment=mol_name + " COMmode " + str(i))

    if debug:
        print("INERTIA: xyz[:5] is", xyz[:5])
    I = inertia_tensor(xyz, mass_per_atom)
    if debug:
        print("INERTIA: xyz[:5] is", xyz[:5])
    w,X = np.linalg.eigh(I)
    if debug:
        print("Interia eigenvals:", w)

    M = np.diag(1.0/(np.repeat(np.sqrt(mass_per_atom), 3)))
    #M = np.diag((np.repeat(np.sqrt(mass_per_atom), 3)))


    # if debug:
    #     print(hess.shape, mass_mat.shape)

    # this is the transformation matrix 
    D = trans_rot_modes(xyz, mass_per_atom, X, rot=True)
    if debug:
        print("D dots:")
        for i in range(D.shape[0]):
            print(np.linalg.norm(D[i]))
    D /= np.linalg.norm(D, axis=1)[:,np.newaxis]
    if debug:
        print("D dots after normalization:")
        for i in range(D.shape[0]):
            print(np.linalg.norm(D[i]))
    #D /= np.linalg.norm(D, axis=0)
    #for i,d in enumerate(D):
    #    norm = np.dot(d,d.T)
    #    if(norm > 1e-7):
    #        D[i] /= np.sqrt(norm)
    N = D.shape[0]
    D = D.T
    TR = D
    if debug:
        print("TR shape")
        print(TR.shape)
        print("D trans/rot dots:")
        print(D)
        
        print("D[:,0] pre qr")
        print(D[:,0])
    if debug:
        testmat = (np.dot(np.dot(TR.T, hess), TR))
        print("TESTMAT TR' * H * TR")
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()
    if debug:
        testmat = np.dot( TR.T, TR)
        print("TR PRE QR DOTS")
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()
    D, _ = np.linalg.qr( D, mode='complete')
    Q, _ = np.linalg.qr( D[:,:N], mode='complete')
    
    # D = np.round(D, 12)
    # Q = np.round(Q, 12)

    for i in range(D.shape[0]):
        if np.linalg.norm( D[:,i]) < 0.0:
            D[:, i] *= -1.0

    if debug:
        print("Q VS D DOT")
        printortho( np.dot( Q.T, D))
        printmat( np.dot( Q.T, D), 2)
    #D /= np.linalg.norm( D, axis=0)
    if debug:
        TR = D[ :, :N]
        testmat = np.dot( TR.T, TR)
        print("TR AFTER QR DOTS")
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()
    if debug:
        testmat = (np.dot(np.dot(TR.T, hess), TR))
        print("TESTMAT")
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()
    #D /= np.linalg.norm(D, axis=1)[:,np.newaxis]
    #D = D.T
    #D = np.array(gram_schmidt(D.T))
    #D = gs2(D.T).T
    if debug:
        print("D shape immediately post qr")
        print(D.shape)
    #D *= -1
    #D = D[6:]
    #Dp = np.eye(D.shape[1])
    #Dp[:6][:] = D
    #D = D.T
    q = np.repeat(np.sqrt(mass_per_atom), xyz.shape[1])
    if debug:
        print("Q (sqrt(m)):")
        print(q)
        print("Dot first N with qr first N")
        for i in range(N):
            x = D[:,i]
            y = TR[:,i]
            print(np.dot(x/np.linalg.norm(x),y/np.linalg.norm(y)))

        print("np.linalg.norm(D,axis=0)")
        print(np.linalg.norm(D,axis=0))
        print("np.linalg.norm(D,axis=1)")
        print(np.linalg.norm(D,axis=1))

        print("D[:,0] post qr")
        print(D[:,0])

        print("D.shape post qr")
        print(D.shape)
        print("full hess is")
        print(hess.shape)

    #if D[:,0].sum() < 0:
    #    D *= -1
    #testmat = np.dot( np.dot( L.T, np.dot( np.dot( D.T, hess / mass_mat), D)), L)
    if debug:
        freq_ic,L = np.linalg.eigh(np.dot(np.dot(D.T, hess), D))
        testmat = np.dot( np.dot( L.T, hess), L)
        print("TEST L D[:10] DOT")
        for row in testmat[:10]:
            for val in row[:10]:
                print("{:11.8f} ".format( val), end="")
            print()
    if debug:
        testmat = np.dot( np.dot( D.T, hess), D)
        print("TEST D H D DOT FULL")
        for row in testmat[:10]:
            for val in row[:10]:
                print("{:11.8f} ".format( val), end="")
            print()
        print("EIGS FULL:", converteig(freq_ic[:10]) )
    if not stdtr:
        _,D = np.linalg.eigh(hess)
    if remove > 0:
        TR = np.array(D[:,:min(N, remove)])
        D = np.array(D[:,min(N, remove):])
        if debug:
            print("TR.shape")
            print(TR.shape)
    freq_ic,L = np.linalg.eigh(np.dot(np.dot(D.T, hess), D))

    testmat = np.dot( np.dot( L.T, np.dot( np.dot( D.T, hess), D)), L)
    a,b = np.linalg.eigh( testmat)
    freq_mol = a

    if debug:
        print("AFTER REMOVE")
        print("TEST L D[:10] DOT")
        for row in testmat[:10]:
            for val in row[:10]:
                print("{:11.8f} ".format( val), end="")
            print()
    if debug:
        testmat = np.dot( np.dot( D.T, hess), D)
        print("EIGS LDHDL:", converteig(a))
        print("TEST D H D DOT FULL")
        for row in testmat[:10]:
            for val in row[:10]:
                print("{:11.8f} ".format( val), end="")
            print()
        print("EIGS FULL:", converteig(freq_ic[:10]))
    
    #q = np.full(xyz.shape, bohr2angstrom) * np.sqrt(mass_per_atom[:,np.newaxis])
    #q = q.reshape(-1)
    if debug:
        print("q")
        print(q)
    S = np.dot(D.T,q)
    if debug:
        print("S[0:21]")
        print(S[0:21])
        print("np.dot(D.T, hess, D) shape")
        print(np.dot(np.dot(D.T, hess), D).shape)

    # Get the modes in IC, which is the transformed Hess
    if debug:
        freq_ic,L = np.linalg.eigh(np.dot(np.dot(D.T, hess), D))
        print("freq_ic")
        print(np.sign(freq_ic)*np.sqrt(np.abs(freq_ic))*conversion)

    
    if debug:
        print(M.shape, D.shape, L.shape)

    mol_modes = np.dot(np.dot(M, D), L)
    DL = np.dot(D, L)

    if debug:
        print("mol_modes.shape")
        print(mol_modes.shape)

    try:
        freq_TR,Ltr = np.linalg.eigh(np.dot(np.dot(TR.T, hess), TR))
    except np.linalg.LinAlgError as ex:
        print(ex)
        if return_DL:
            return None, None, None
        else:
            return None, None

    #testmat = np.dot( Ltr, TR.T)
    testmat = np.dot( np.dot( Ltr.T, np.dot( np.dot( TR.T, hess), TR)), Ltr)
    if debug:
        print("TEST Ltr' TR' H TR Ltr DOT (shape):", Ltr.shape, TR.shape)
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( converteig(val) ), end="")
            print()
    a,b = np.linalg.eigh( testmat)
    testmat = b
    freq_TR = a

    if debug:
        print("EIGS TESTMAT:", converteig(a))
        print("EIGS of EIGS")
        for row in testmat:
            for val in row:
                print("{:11.8f} ".format( val), end="")
            print()

    if debug:
        print("freq_TR")
        print(np.sign(freq_TR)*np.sqrt(np.abs(freq_TR))*conversion)

    
        print(M.shape, TR.shape, L.shape)

    TR_modes = np.dot(np.dot(M, TR), Ltr)

    if debug:
        print("TR_modes.shape")
        print(TR_modes.shape)



    if remove > 0:
        mol_modes = np.hstack((TR_modes, mol_modes))
        freq_ic = np.hstack((freq_TR, freq_mol))
        mol_modes[:,:min(remove, N)] = TR_modes
    if debug:
        print("EIGS FULL AFTER COMBINE:", converteig(freq_ic[:10]))
    #else:
    #    Norm = 1./np.linalg.norm(TR, axis=0)

    #    print(TR_modes.shape, Norm.shape)
    #    for i,_ in enumerate(Norm):
    #        TR_modes[:,i] *= Norm[i]


    #for i in range( TR_modes.shape[1]):
    #    TR_modes[:,i] *=  q

    #Norm = np.linalg.norm(mol_modes, axis=0)

    #print(mol_modes.shape, Norm.shape)
    #for i,_ in enumerate(Norm):
    #    mol_modes[:,i] /= Norm[i]

    freq_ic = converteig(freq_ic)
    if verbose:
        for i, mode in enumerate(mol_modes.T):
            with open(mol_name+".mode_"+str(i)+".pure.xyz", 'w') as fid:
                mode_to_xyz(syms, xyz, (mode).reshape(-1, 3) , outfd=fid, 
                    comment=mol_name + " mode " + str(i)  + f" {freq_ic[i]} cm-1")

    if debug:
        raw_mode = mol_modes.copy()
        for i in range( raw_mode.shape[1]):
            raw_mode[:,i] = raw_mode[:,i] / q
        print("cart modes norm (shape=", mol_modes.shape, ":")
        print(np.linalg.norm(mol_modes, axis=0))
        print("cart modes orthogonal (eps=1e-8)? :")
        orth = np.dot( mol_modes.T, mol_modes)
        for i in orth:
            for x in i:
                x = abs(x)
                if x > 1e-8:
                    print("{:3.1f}".format(x), end="")
                else:
                    print("Y", end="")
            print()

        print("Angular momentum of cart:")
        for i, mode in enumerate(mol_modes.T):
            print( "mode {:3d}: {:4.2e}".format( i, 
                angular_momentum( xyz, mass_per_atom, mode)))

        print("raw modes pre-norm:")
        print(np.linalg.norm( raw_mode, axis=0))
        print("raw modes orthogonal (eps=1e-8)? :")
        raw_mode /= np.linalg.norm( raw_mode, axis=0)
        print("raw modes post-norm:")
        print(np.linalg.norm( raw_mode, axis=0))
        orth = np.dot( raw_mode.T, raw_mode)
        for i in orth:
            for x in i:
                x = abs(x)
                if x > 1e-8:
                    print("{:3.1f}".format(x), end="")
                else:
                    print("Y", end="")
            print()

        print("Angular momentum of raw:")
        for i, mode in enumerate(raw_mode.T):
            print( "mode {:3d}: {:4.2e}".format( i, 
                angular_momentum( xyz, mass_per_atom, mode)))
            if verbose:
                with open(mol_name+".mode_"+str(i)+".raw.xyz", 'w') as fid:
                    mode_to_xyz(syms, xyz, (mode).reshape(-1, 3) , outfd=fid, 
                            comment=mol_name + " mode " + str(i))

    if return_DL:
        return freq_ic, mol_modes, DL
    else:
        return freq_ic, mol_modes
