from bson.objectid import ObjectId
import math
import numpy as np
from scipy.stats import norm
import sys


N_CODE = 0
A_CODE = 1
B_CODE = 2
H_CODE = 3


def _solve_quad(a, b, c):
    """
    Find solutions for the given quadratic equation where
    :param a: quadratic coefficient
    :param b: linear coefficient
    :param c: constant term
    :return: the pair of values for x that solves ax^2 + bx + c = 0
    """
    # TODO there are more numerically stable implementations you should look into using
    sqrt_term = np.sqrt(b ** 2 - 4 * a * c)
    x1 = (-b + sqrt_term) / (2 * a)
    x2 = (-b - sqrt_term) / (2 * a)

    return x1, x2


def prob_density(cluster, point_x, point_y):
    """
    Find the probability density of the given cluster at x, y
    :param cluster: the cluster
    :param point_x: the x point
    :param point_y: the y point
    :return: the value of the probability density
    """

    mean_x = cluster['mean_x']
    mean_y = cluster['mean_y']
    rot_x_var = cluster['rot_x_var']
    rot_y_var = cluster['rot_y_var']

    clust_rot_mat = inverse_rot(mean_x, mean_y)

    rot_point_x, rot_point_y = clust_rot_mat.dot([point_x, point_y])

    # y is already centered around zero because of the rotation
    # but we still need to center x
    clust_hyp, _ = clust_rot_mat.dot([mean_x, mean_y])
    prob_density_x = norm.pdf(rot_point_x - clust_hyp, scale=np.sqrt(rot_x_var))
    prob_density_y = norm.pdf(rot_point_y, scale=np.sqrt(rot_y_var))

    # because we rotated for diagonal covariance we can simply multiply the
    # x and y covariances
    return prob_density_x * prob_density_y


def inverse_rot(x, y):
    """
    Returns a 2D rotation matrix which will invert thetas rotation
    :param theta: the rotation angle to invert
    :return: the rotation matrix which will inverse thetas rotation
    """

    # a 2D rotation matrix looks like:
    #
    # R = | cos(theta), -sin(theta) |
    #     | sin(theta), cos(theta)  |
    #
    # in order to invert the rotation we can just transpose the matrix giving
    #
    # R = |  cos(theta), sin(theta) |
    #     | -sin(theta), cos(theta) |

    hyp = math.sqrt(x * x + y * y)
    if hyp == 0:
        return np.array([
            [0.0, 0.0],
            [0.0, 0.0],
        ])
    else:
        cos_theta = x / hyp
        sin_theta = y / hyp
        return np.array([
            [cos_theta,  sin_theta],
            [-sin_theta, cos_theta],
        ])


def _scale_vector_to_one(v):
    """
    Scales the given vector in place such that sum(v) == 1
    """
    v /= np.sum(v)


def _scale_matrix_rows_to_one(m):
    """
    Scales the given vector in place such that all np.sum(m, 1) == 1
    """
    # get a column vector of row sums
    row_sums = np.sum(m, 1)
    row_sums = np.reshape(row_sums, (-1, 1))
    m /= row_sums


def samples_to_ab_codes(samples, chromosome, snps):
    """
    Converts the "GACT-" calls into "ABHN" which is a more natural representation for the HMM to work with
    :param samples: the samples from mongo DB whose SNP calls we're going to extract and convert
    :param chromosome: the chromosome to convert
    :param snps: the SNP annotation dicts from mongo DB
    :return: the AB codes organized in as an NP array
    """
    snp_count = 0
    x_calls = []
    y_calls = []
    for snp in snps:
        snp_count += 1
        x_calls.append(snp['x_probe_call'])
        y_calls.append(snp['y_probe_call'])
    x_calls = np.array(x_calls)
    y_calls = np.array(y_calls)
    ab_codes = np.empty((snp_count, len(samples)), dtype=np.uint8)
    ab_codes.fill(255)

    for i, curr_sample in enumerate(samples):
        if 'allele1_fwds' in curr_sample['chromosome_data'][chromosome]:
            allele1_fwds = np.array(curr_sample['chromosome_data'][chromosome]['allele1_fwds'])
            allele2_fwds = np.array(curr_sample['chromosome_data'][chromosome]['allele2_fwds'])

            # here we convert nucleotides (GACT and '-' for no call) into AB codes
            allele1_is_a = allele1_fwds == x_calls
            allele2_is_a = allele2_fwds == x_calls
            allele1_is_b = allele1_fwds == y_calls
            allele2_is_b = allele2_fwds == y_calls
            alleles_are_het = np.logical_or(
                np.logical_and(allele1_is_a, allele2_is_b),
                np.logical_and(allele1_is_b, allele2_is_a))
            ab_codes[np.logical_and(allele1_is_a, allele2_is_a), i] = A_CODE
            ab_codes[np.logical_and(allele2_is_b, allele2_is_b), i] = B_CODE
            ab_codes[alleles_are_het, i] = H_CODE
            ab_codes[np.logical_or(allele1_fwds == '-', allele2_fwds == '-'), i] = N_CODE
        elif 'snps' in curr_sample['chromosome_data'][chromosome]:
            sample_snps = np.array(curr_sample['chromosome_data'][chromosome]['snps'])

            ab_codes[sample_snps == x_calls, i] = A_CODE
            ab_codes[sample_snps == y_calls, i] = B_CODE
            ab_codes[sample_snps == 'H', i] = H_CODE
            ab_codes[sample_snps == '-', i] = N_CODE
        else:
            raise Exception('sample does not have allele1_fwds or snps attributes for chromosome {}'.format(chromosome))

        bad_codes = ab_codes[:, i] == 255
        num_bad_codes = np.count_nonzero(bad_codes)
        if num_bad_codes:
            err_msg = 'found unexpected {} SNP codes in sample: {}, chr: {}'.format(
                    num_bad_codes,
                    curr_sample['sample_id'],
                    chromosome)
            print(err_msg)
            ab_codes[bad_codes, i] = N_CODE
            #raise Exception(err_msg)

    return ab_codes


class SnpHaploHMM:

    def __init__(self, trans_prob, hom_obs_probs, het_obs_probs, n_obs_probs):
        """
        construct a SNP Haplotype HMM. Note that the observation probabilities have a different structure
        than you would see in a standard HMM. This is because the valid observations change on
        a per SNP basis. We define N to be the number of haplotypes
        :param trans_prob:
            the probability of transitioning from one haplotype to another in the space of
            a single SNP
        :param hom_obs_probs:
            a categorical distribution describing P(obs | homozygous read from hidden state).
            The values in this array will be scaled such that sum(hom_obs_probs) == 1.
            hom_obs_probs should be a numerical array where the indexes mean:
            0) observation is homozygous matching the hidden state (either A->A or B->B). For a
               well behaved HMM this should be by far the most probable. The other observations
               would be considered genotyping errors
            1) observation is homozygous not matching the hidden state (either A->B or B->A)
            2) observation is heterozygous (either A->H or B->H)
            3) observation is no read (either A->N or B->N)
        :param het_obs_probs:
            a categorical distribution describing P(obs | heterozygous read from hidden state).
            The values in this array will be scaled such that sum(het_obs_probs) == 1.
            het_obs_probs should be a numerical array where the indexes mean:
            0) observation is heterozygous (H->H). For a well behaved HMM this should be by far
               the most probable. The other observations would be considered genotyping errors
            1) observation is homozygous (either H->A or H->B)
            2) obs is no read (H->N)
        :param n_obs_probs:
            a categorical distribution describing P(obs | no read from hidden state). We'll call
            any of NN, NH, NA or NB an 'N' hidden state (that is to say that only one of the two
            contributing genotypes needs to be N for us to call the hidden state N)
            The values in this array will be scaled such that sum(n_obs_probs) == 1.
            n_obs_probs should be a numerical array where the indexes mean:
            0) observation is no read (N->N)
            1) observation is homozygous (N->A or N-B)
            2) observation is heterozygous (N->H)
        """

        self.trans_prob = trans_prob

        # scale and assign observation probs
        _scale_vector_to_one(hom_obs_probs)
        self.hom_obs_probs = hom_obs_probs
        _scale_vector_to_one(het_obs_probs)
        self.het_obs_probs = het_obs_probs
        _scale_vector_to_one(n_obs_probs)
        self.n_obs_probs = n_obs_probs

        # pack the sparse observation probs into a dense matrix for quick lookup
        # using codes where row is hidden state encoded as 0->N, 1->A, 2->B, 3->H and column
        # is observation encoded as 0->N, 1->A, 2->B, 3->H
        self.obs_prob_matrix = np.zeros((4, 4, 4), dtype=np.float64)
        for hap1_call in range(4):
            for hap2_call in range(4):
                for obs_call in range(4):
                    # if no read from hidden state (we consider it a "no read" if either
                    # parent is N or H)
                    if hap1_call == N_CODE or hap2_call == N_CODE or hap1_call == H_CODE or hap2_call == H_CODE:
                        # if N from observation
                        if obs_call == N_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.n_obs_probs[0]
                        # if hom from observation
                        elif obs_call == A_CODE or obs_call == B_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.n_obs_probs[1]
                        # else het from observation
                        else:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.n_obs_probs[2]
                    # if hom from hidden state
                    elif hap1_call == hap2_call:
                        # if N from observation
                        if obs_call == N_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.hom_obs_probs[3]
                        # if het from observation
                        elif obs_call == H_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.hom_obs_probs[2]
                        # if matching hom from observation
                        elif obs_call == hap1_call:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.hom_obs_probs[0]
                        # else non-matching hom from observation
                        else:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.hom_obs_probs[1]
                    # else het from hidden state
                    else:
                        # if N from observation
                        if obs_call == N_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.het_obs_probs[2]
                        # if het from observation
                        elif obs_call == H_CODE:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.het_obs_probs[0]
                        # else hom from observation (A or B)
                        else:
                            self.obs_prob_matrix[hap1_call, hap2_call, obs_call] = self.het_obs_probs[1]

    def viterbi(self, haplotype_ab_codes, observation_ab_codes):
        """
        Run the viterbi algorithm to find the haplotype state sequence (S) that maximizes
        P(S|O, HMM) where O is the observation sequence and HMM is this parameterized HMM
        :param haplotype_ab_codes:
            a 2D matrix where row_index == snp_index and col_index == haplotype_index (which must
            match the haplotype indices used in the transition probs).
            The numerical genotype codes used should be: 0->N, 1->A, 2->B, 3->H
        :param observation_ab_codes:
            a vector of genotype observations the length of observation_ab_codes should match
            the row count of haplotype_ab_codes
            The numerical genotype codes used should be: 0->N, 1->A, 2->B, 3->H
        :return:
            the tuple (max_likelihood_states, log_likelihood) where,
            * max_likelihood_states is the most likely sequence of states given the HMM and the observations.
              This will be a numeric vector where indices match observation_ab_codes and values represent
              haplotype indices
            * log_likelihood is the log likelihood of this hidden state sequence
        """
        obs_count, haplotype_count = haplotype_ab_codes.shape

        state_hap1_indices, state_hap2_indices = np.triu_indices(haplotype_count)
        state_count = state_hap1_indices.size

        # initialize transition probabilities
        log_transition_to_new_state_prob = 0
        if state_count > 1:
            log_transition_to_new_state_prob = np.log(self.trans_prob / (state_count - 1))

        log_transition_to_same_state_prob = np.log(1.0 - self.trans_prob)
        log_trans_probs = np.empty((state_count, state_count), dtype=np.float64)
        log_trans_probs[:] = log_transition_to_new_state_prob
        np.fill_diagonal(log_trans_probs, log_transition_to_same_state_prob)

        # set uniform initial probs
        log_init_prob = np.log(1.0 / state_count)
        log_obs_prob_matrix = np.log(self.obs_prob_matrix)

        # initialize viterbi using first observations and the init_probs
        curr_log_obs_probs = log_obs_prob_matrix[
            haplotype_ab_codes[0, state_hap1_indices],
            haplotype_ab_codes[0, state_hap2_indices],
            observation_ab_codes[0]]
        curr_log_likelihoods = log_init_prob + curr_log_obs_probs

        # step forward in time updating the "from" state likelihoods for each
        # "to" state at each step
        from_state_lattice = np.zeros((obs_count - 1, state_count), dtype=np.uint16)
        for t in range(1, obs_count):
            prev_log_likelihoods = curr_log_likelihoods
            curr_log_likelihoods = np.zeros(state_count, dtype=np.float64)

            for to_state in range(state_count):
                # TODO is log_trans_probs[:, to_state] right??
                from_log_likelihoods = prev_log_likelihoods + log_trans_probs[:, to_state]
                max_from_state = np.argmax(from_log_likelihoods)
                from_state_lattice[t - 1, to_state] = max_from_state
                curr_log_likelihoods[to_state] = from_log_likelihoods[max_from_state]

            curr_log_obs_probs = log_obs_prob_matrix[
                haplotype_ab_codes[t, state_hap1_indices],
                haplotype_ab_codes[t, state_hap2_indices],
                observation_ab_codes[t]]
            curr_log_likelihoods += curr_log_obs_probs

        # backtrace through the most likely path starting with the final state
        max_final_state = np.argmax(curr_log_likelihoods)
        max_final_likelihood = curr_log_likelihoods[max_final_state]
        max_likelihood_states = np.zeros(obs_count, dtype=np.uint16)
        max_likelihood_states[obs_count - 1] = max_final_state
        for t in reversed(range(obs_count - 1)):
            max_likelihood_states[t] = from_state_lattice[t, max_likelihood_states[t + 1]]

        max_likelihood_states = [
            (int(state_hap1_indices[s]), int(state_hap2_indices[s]))
            for s in max_likelihood_states
        ]
        return max_likelihood_states, max_final_likelihood

    def log_likelihood(self, haplotype1_ab_codes, haplotype2_ab_codes, observation_ab_codes):
        """
        Compute the log likelihood of the observation_ab_codes sequece of SNPs given that the parental
        strains' genotypes are represented by haplotype1_ab_codes and haplotype2_ab_codes
        :param haplotype1_ab_codes: sequence of SNPs from the first parental haplotype converted to AB codes
        :param haplotype2_ab_codes: sequence of SNPs from the second parental haplotype converted to AB codes
        :param observation_ab_codes: sequence of SNPs from the sample we're haplotyping
        :return: the log probability
        """
        log_obs_prob_matrix = np.log(self.obs_prob_matrix)
        return np.sum(log_obs_prob_matrix[haplotype1_ab_codes, haplotype2_ab_codes, observation_ab_codes])

    def forward_scan(self, haplotype_ab_codes, observation_ab_codes):
        """
        :param haplotype_ab_codes:
        :param observation_ab_codes:
        :return:
        """

        # forward_likelihoods = None
        # forward_scaling_factors = None
        #
        # return forward_likelihoods, forward_scaling_factors
        raise Exception('implement me')


# def main():
#
#     ###########################################################################
#     # This main function is just some smoke-signal test code for this module. #
#     ###########################################################################
#
#     hom_obs_probs = np.array([
#         50,     # matching homozygous
#         0.5,    # opposite homozygous
#         1,      # het
#         1,      # no read
#     ], dtype=np.float64)
#     het_obs_probs = np.array([
#         50,     # het obs
#         2,      # hom obs
#         2,      # no read obs
#     ], dtype=np.float64)
#     n_obs_probs = np.ones(3, dtype=np.float64)
#
#     trans_prob = 0.001
#     test_hmm = SnpHaploHMM(trans_prob, hom_obs_probs, het_obs_probs, n_obs_probs)
#
#     # Each sample is a different set of parental strains but only 129S, B6J, FVB for the most part
#     print('getting AB Codes')
#     for chr in (str(i) for i in range(1, 20)):
#         print('chromosome', chr)
#         ab_codes = sample_ids_to_ab_codes(
#             [
#                 # only 'C57BL/6NJm39418', '129S7m', 'FVB001', are expected to really be parental strains
#                 'UM-004', 'C57BL/6NJm39418', '129S7m', 'FVB001', 'C3H/HeNTac', 'BALB/cByJm',
#                 'NOR/LtJm32299', 'NON/ShiLtJm38111', 'SJL/Jm35807', 'ALR/LtJm34133', 'CZECHII/EiJ',
#                 'ZALENDE/EiJ', 'SPRET/EiJ', 'ALR/LtJm34133',
#             ],
#             chr,
#         )
#         print('performing viterbi')
#         np.set_printoptions(threshold=np.nan)
#         max_likelihood_states, max_final_likelihood = test_hmm.viterbi(
#             haplotype_ab_codes=ab_codes[:, 1:],
#             observation_ab_codes=ab_codes[:, 0])
#         print(max_likelihood_states)
#         print(max_final_likelihood)
#         print('done with viterbi')
#
# if __name__ == '__main__':
#     main()
#
##################### Robyn's code
#
#
import pandas

haplotype_ab_codes = pandas.read_csv('haplotype_ab_codes.csv').to_numpy()
observation_ab_codes = pandas.read_csv('observation_ab_codes.csv').to_numpy()

hom_obs_probs = np.array([0.97, 0.01,0.01, 0.01], dtype=np.float64)
het_obs_probs = np.array([0.98,0.01,0.01], dtype=np.float64)
trans_prob = 0.001
n_obs_probs = np.ones(3, dtype=np.float64)

test_hmm = SnpHaploHMM(trans_prob, hom_obs_probs, het_obs_probs, n_obs_probs)
max_likelihood_states, max_final_likelihood = test_hmm.viterbi(haplotype_ab_codes = haplotype_ab_codes, observation_ab_codes = observation_ab_codes)

np.savetxt("max_likelihood_states.csv", max_likelihood_states, delimiter=",")
