"""
This script checks that PMBEC-as-a-distance-measure is a metric space.

As it turns out, it is a metric space if we set MAX = maximum value of
the PMBEC convariance matrix, but it is not a metric space if we set
MAX = some local maximum (e.g. the max of the row and column).

Given that residue distance is metric, it follows that peptide distance
is metric: but this script also checks that (using randomly generated
3-mers and randomly generated positional weights).
"""

import random
from immuno.pmbec.distance import PMBECDistance

DEBUG = False

num_failures = 0
distance = PMBECDistance([1] * 9)
for key_a in distance.covariance.keys():
    for key_b in distance.covariance.keys():
        row_a = key_a[0]
        row_b = key_b[0]
        col_a = key_a[1]
        col_b = key_b[1]
        # Example: key_a = AV, key_b = VK, transitive = dist(A, K)
        if col_a == row_b:
            transitive_distance = distance._residue_distance(row_a, col_b)
            a_distance = distance._residue_distance(row_a, col_a)
            b_distance = distance._residue_distance(row_b, col_b)
            added_distance = a_distance + b_distance
            if transitive_distance > added_distance:
                num_failures += 1
                print("Triangle failure: " + key_a + ", " + key_b)
            else:
                print("Success: %s (%f) + %s (%f) = %f >= %s%s (%f)" %
                    (key_a, a_distance, key_b, b_distance, added_distance,
                     row_a, col_b, transitive_distance))
print("Number of triangle failures: %d" % num_failures)

def get_random_peptide():
    peptide_string = ""
    for i in xrange(3):
        peptide_string += random.choice(AMINO_ACIDS)
    return peptide_string

AMINO_ACIDS = ['W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
    'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q']
num_failures = 0
for i in xrange(0, 1000000):
    distance = PMBECDistance([random.random() for _ in xrange(3)],
            peptide_length = 3)
    peptide_string_a = get_random_peptide()
    peptide_string_b = get_random_peptide()
    peptide_string_c = get_random_peptide()
    distance_ab = distance.get_distance(peptide_string_a, peptide_string_b)
    distance_bc = distance.get_distance(peptide_string_b, peptide_string_c)
    distance_ac = distance.get_distance(peptide_string_a, peptide_string_c)
    added_distance = distance_ab + distance_bc
    if distance_ac - added_distance >= 0.0000000001:
        num_failures += 1
        print("Triangle failure: %s, %s, %s / %f, %f, %f, %f" %
                (peptide_string_a, peptide_string_b, peptide_string_c,
                 distance_ab, distance_bc, added_distance, distance_ac))
    else:
        if (DEBUG):
            print("Success: %s, %s, %s / %f, %f, %f, %f" %
                    (peptide_string_a, peptide_string_b, peptide_string_c,
                     distance_ab, distance_bc, added_distance, distance_ac))
print("Number of peptide triangle failures: %d" % num_failures)
