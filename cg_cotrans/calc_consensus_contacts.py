# Copyright (C) 2017 William M. Jacobs

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
import sys

import argparse, prody, gzip, subprocess
import numpy as np
from itertools import product
from Bio import pairwise2

from codons import codon_to_aa
from polymer import Polymer
from substructures import find_substructures, calc_co

def aaseq(seq):
    codons = [seq[3*i:3*(i+1)] for i in range(len(seq)//3)]
    return ''.join(codon_to_aa[c] for c in codons if codon_to_aa[c] != 'Stop')

def align_sequences(seqs, clustalo_cmd='/usr/bin/clustalo'):
    with subprocess.Popen([clustalo_cmd, '-i', '-'],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) \
                          as proc:
        clustalo_input = bytes(''.join(">%d\n%s\n" % (i, seqs[i]) for i in range(len(seqs))), \
                               'ascii')
        clustalo_output, clustalo_err = proc.communicate(clustalo_input)
    aligned_seqs = []
    for line in clustalo_output.decode('ascii').split('\n'):
        if len(line) > 0 and line[0] == '>':
            aligned_seqs.append('')
        elif len(line) > 0:
            aligned_seqs[-1] += line.strip()
    return aligned_seqs

def get_alignment_indices(aligned_seqs):
    indices = []
    for i in range(1, len(aligned_seqs)):
        indices.append({})
        c = 0
        for j in range(len(aligned_seqs[i])):
            if aligned_seqs[i][j] != '-':
                indices[-1][c] = j
                c += 1
    return indices

def get_percent_id(aligned_seqs):
    percid = []
    for i in range(1, len(aligned_seqs)):
        percid.append(sum(1 if aligned_seqs[0][c] == aligned_seqs[i][c] else 0 \
                          for c in range(len(aligned_seqs[0])) if aligned_seqs[0][c] != '-') \
                      / len(aligned_seqs[0].replace('-', '')) * 100.)
    return percid

def get_shifted_residues(chain):
    residues = [r for r in chain.iterResidues() if r.isprotein]
    shifted_residues = [[] for i in range(len(residues))]
    for i in range(len(residues)):
        residue = residues[i]
        for atom in residue.iterAtoms():
            if atom.getElement() == 'H':
                continue
            elif atom.getName() == 'N' and i > 0:
                shifted_residues[i - 1].append(atom)
            else:
                shifted_residues[i].append(atom)
    return [res for res in shifted_residues if len(res) > 0]

def build_contact_map_from_residues(n, indices, residues, cutoff=4., min_sequence_distance=3):
    atoms = np.array([[x for x in atom.getCoords()] for residue in residues \
                      for atom in residue if atom.getElement() != 'H'])
    atomranges = []
    j = 0
    for residue in residues:
        start = j
        for atom in residue:
            if atom.getElement() != 'H':
                j += 1
        stop = j
        atomranges.append((start, stop))
    centers = np.array([np.mean([atoms[ai] for ai in range(*atomranges[i])]) \
                        for i in range(len(residues))])
    A = np.zeros((n, n))
    for i in range(len(residues)):
        for j in range(i + min_sequence_distance, len(residues)):
            if np.linalg.norm(centers[i] - centers[j]) <= cutoff * 4:
                for ai in range(*atomranges[i]):
                    for aj in range(*atomranges[j]):
                        if np.linalg.norm(atoms[ai] - atoms[aj]) <= cutoff:
                            A[indices[i],indices[j]] += 1
                            A[indices[j],indices[i]] += 1
    return A

def identify_backbone_hydrogen_bonds(contacts, indices, residues, cutoff=4.):
    H = np.zeros((len(contacts), len(contacts)))
    for i in indices:
        for j in indices:
            if contacts[indices[i],indices[j]] > 0:
                for atom1, atom2 in product(residues[i], residues[j]):
                    if ((atom1.getName() == 'N' and atom2.getName() == 'O') or \
                        (atom1.getName() == 'O' and atom2.getName() == 'N')) and \
                        np.linalg.norm(atom1.getCoords() - atom2.getCoords()) <= cutoff:
                        H[indices[i],indices[j]] = H[indices[j],indices[i]] = 1
                        break
    return H

def consensus_contacts(wtseq, allcontacts, allhbonds, consensus_frac=0.25):
    n = list(allcontacts.values())[0].shape[0]
    wtn = len(wtseq.replace('-', ''))
    contacts = np.zeros((wtn, wtn))
    hbonds = np.zeros((wtn, wtn))
    ii = -1
    for i in range(n):
        if wtseq[i] != '-':
            ii += 1
            jj = -1
            for j in range(n):
                if wtseq[j] != '-':
                    jj += 1
                    Cij = np.array([allcontacts[pdbid][i,j] for pdbid in allcontacts.keys()])
                    Hij = np.array([allhbonds[pdbid][i,j] for pdbid in allcontacts.keys()])
                    if sum(1 if Cij[k] > 0 else 0 for k in range(len(Cij))) / len(Cij) \
                       >= consensus_frac:
                        contacts[ii,jj] = np.max(Cij)
                        hbonds[ii,jj] = np.max(Hij)
    return contacts, hbonds

def bond_energies(polymer, substructures, contacts, hbonds, alpha_helix=0.625, alpha_hbond=16.):
    ishelix = {}
    for i in range(contacts.shape[0]):
        for j in range(contacts.shape[1]):
            ishelix[(i,j)] = ishelix[(j,i)] = 0
            for ss in substructures:
                if (i,j) in ss or (j,i) in ss:
                    if calc_co(ss) < 5:
                        ishelix[(i,j)] = ishelix[(j,i)] = 1
                        break
    U = np.zeros(contacts.shape)
    for i in range(contacts.shape[0]):
        for j in range(contacts.shape[1]):
            U[i,j] = alpha_helix**ishelix[(i,j)] * (contacts[i,j] / alpha_hbond + hbonds[i,j])
    return U

def write_polymer(U, stream):
    stream.write("nresidues = %d\n" % U.shape[0])
    for i in range(U.shape[0]):
        for j in range(i, U.shape[1]):
            if U[i,j] > 0:
                stream.write("%d %d %g\n" % (i, j, U[i,j]))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, help="path to fasta sequence")
    parser.add_argument('--pdbids', type=str, help="comma-separated list of pdbids")
    parser.add_argument('--output', type=str, metavar='PATH', default='polymer.dat', \
                        help="path to output file [polymer.dat]")
    parser.add_argument('--structure-path', type=str, metavar='PATH', default='./', \
                        help="path to gzip'ed input pdb files [./]")
    parser.add_argument('--min-percent-id', type=float, default=25., \
                        help="minimum percent sequence identity [95.]")
    parser.add_argument('--consensus-frac', type=float, default=0.25, \
                        help="minimum contact consensus fraction [0.25]")
    clargs = parser.parse_args()
###
    #print("calc_consensus...")
    #sys.exit()
###


    seq = ''
    with open(clargs.fasta, 'r') as f:
        for line in f:
            if '>' in line:
                continue
            else:
                seq += line.strip()
    wtseq = aaseq(seq)

    chains = {}
    for pdbid in clargs.pdbids.split(','):
        pdbfile = clargs.structure_path + '%s.pdb.gz' % pdbid
        with gzip.open(pdbfile, 'rt') as f:
            try:
                mol, header = prody.parsePDBStream(f, header=True)
            except ValueError as e:
                print("PDB Error: %s: %s" % (pdbid, e))
                continue
        if 'experiment' not in header or 'X-RAY' not in header['experiment']:
            print("%s is not an X-ray structure; skipping" % clargs.pdbfile)
        else:
            if 'resolution' in header:
                resolution = float(header['resolution'])
            else:
                resolution = 100.
            pdbchains = {}
            for chain in mol.iterChains():
                chainseq = chain.getSequence()
                score = pairwise2.align.globalms(wtseq, chainseq, 2, -1, -1, -.5, score_only=True)
                if isinstance(score, list):
                    continue
                residues = get_shifted_residues(chain)
                if len(residues) != len(chainseq):
                    raise Exception("mismatched lengths")
                pdbchains[str(chain).split()[-1]] = score, (chainseq, residues, score)
            bestchain = max(pdbchains, key=lambda y: pdbchains[y][0])
            chains[pdbid] = pdbchains[bestchain][1]
    if len(chains) > 100:
        chains = {pdbid : chains[pdbid] for pdbid in sorted(chains, key=lambda y: -chains[y][-1])[:100]}

    pdbids = sorted(chains.keys())
    if len(pdbids) == 0:
        print("No chains found.")
        raise SystemExit
    print("Found %d chains." % len(chains))
    aligned_chains = align_sequences([wtseq] + [chains[pdbid][0] for pdbid in pdbids])
    aligned_indices = get_alignment_indices(aligned_chains)
    n = len(aligned_chains[0])

    contact_maps = {}
    hbonds = {}
    for i in range(len(pdbids)):
        indices = aligned_indices[i]
        residues = chains[pdbids[i]][1]
        contact_maps[pdbids[i]] = build_contact_map_from_residues(n, indices, residues)
        hbonds[pdbids[i]] = identify_backbone_hydrogen_bonds(contact_maps[pdbids[i]], indices, residues)

    percid = get_percent_id(aligned_chains)
    for i in range(len(pdbids)):
        print("%s %6.2f" % (pdbids[i], percid[i]))
    selected_pdbids = [pdbids[i] for i in range(len(pdbids)) if percid[i] >= clargs.min_percent_id]
    if len(selected_pdbids) == 0:
        print("No chains selected.")
        raise SystemExit
    print("Using %d chains with percentid >= %g." % (len(selected_pdbids), clargs.min_percent_id))

    contacts, hbonds = consensus_contacts(aligned_chains[0], \
                                          {pdbid : contact_maps[pdbid] for pdbid in selected_pdbids}, \
                                          {pdbid : hbonds[pdbid] for pdbid in selected_pdbids}, \
                                          consensus_frac=clargs.consensus_frac)
    polymer = Polymer([(i,j) for i in range(contacts.shape[0]) for j in range(contacts.shape[1]) \
                       if contacts[i,j] > 0])
    substructures = find_substructures(polymer)
    U = bond_energies(polymer, substructures, contacts, hbonds)
    print("Number of residues:", contacts.shape[0])
    print("Minimum consensus fraction:", clargs.consensus_frac)
    print("Total bond energy:", 0.5 * U.sum())

    print("Writing %s" % clargs.output)
    with open(clargs.output, 'w') as f:
        write_polymer(U, f)
