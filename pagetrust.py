import numpy as np


def pagetrust(
    G,
    alpha=0.85,
    M=0,
    beta=1.0,
    z=None,
    max_iter=100,
    tol=1e-6,
    forced_zero_nodes=None,
    handle_dangling=True,
):
    """
    Compute the PageTrust scores for a directed graph G.

    Parameters:
    - G: A directed graph (networkx.DiGraph) where edges can have positive or negative weights.
    - alpha: Damping factor for the PageTrust algorithm (default: 0.85).
    - M: A parameter that controls the influence of the teleportation vector (default: 0).
    - beta: A parameter that controls the influence of the negative edges (default: 1.0).
    - z: Teleportation vector (default: uniform distribution).
    - max_iter: Maximum number of iterations for convergence (default: 100).
    - tol: Tolerance for convergence (default: 1e-6).
    - forced_zero_nodes: A list of nodes for which the trust score should be forced to zero (default: None).
    - handle_dangling: If True, handle dangling nodes by redistributing their trust mass according to z (default: True).
    """

    # Create a mapping from node to index for matrix representation
    nodes = list(G.nodes())
    n = len(nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    # Create a mask for nodes that should have their trust score forced to zero
    forced_zero_mask = np.zeros(n, dtype=bool)
    if forced_zero_nodes:
        for node in forced_zero_nodes:
            if node in node_to_idx:
                forced_zero_mask[node_to_idx[node]] = True

    # Create adjacency matrices for positive and negative edges
    A_pos = np.zeros((n, n))
    A_neg = np.zeros((n, n))

    # Fill the adjacency matrices based on edge weights
    for u, v, data in G.edges(data=True):
        weight = data.get("weight", 1.0)
        i = node_to_idx[u]
        j = node_to_idx[v]
        if weight > 0:
            A_pos[i, j] = weight
        elif weight < 0:
            A_neg[i, j] = 1.0  # L^- indicator

    # Normalize the teleportation vector z
    if z is None:
        z = np.ones(n) / n
    else:
        z = np.array(z, dtype=float)
        z /= np.sum(z)

    # If there are nodes that should have their trust score forced to zero, we set their corresponding entries in z to zero and renormalize z.
    if forced_zero_mask.any():
        z = z.copy()
        z[forced_zero_mask] = 0.0
        z_sum = z.sum()
        if z_sum > 0:
            z /= z_sum

    # This is a new addition beyond the original paper.
    # The original PageTrust paper does not handle dangling nodes
    # (nodes with no positive out-links)
    # To avoid this problem we replace the row of A_pos for dangling nodes with
    # z, which ensures that the trust mass leaving a dangling node is
    # distributed according to z, rather than being trapped in a self-loop.
    # Effectively, this changes how the node passes it's trust value to other nodes,
    # but has no bearing on the trust of the node itself.
    out_degrees = np.sum(A_pos, axis=1)
    dangling_mask = out_degrees == 0
    A_pos_eff = A_pos.copy()
    if dangling_mask.any() and handle_dangling:
        A_pos_eff[dangling_mask, :] = z
    out_degrees_eff = np.sum(A_pos_eff, axis=1)

    # Initialize the trust scores and the P matrix
    x = np.ones(n) / n
    P = np.copy(A_neg)
    P_tilde_diag = np.diag(P).copy()

    # Iterate until convergence or until the maximum number of iterations is reached
    for iteration in range(max_iter):
        print(f"Iteration {iteration + 1}")
        # Store the previous trust scores for convergence checking
        x_prev = np.copy(x)
        # Compute the diagonal matrix D_P_beta based on the current P_tilde_diag and beta
        D_P_beta = np.power(np.clip(1.0 - P_tilde_diag, 0, 1), beta)
        # Compute the transition matrix T based on the current trust scores and the adjacency matrices
        T = np.zeros((n, n))
        # Compute the incoming positive sum for each node, which is used in the denominator of the transition matrix calculation
        incoming_pos_sum = np.zeros(n)
        for i in range(n):
            for k in range(n):
                if A_pos_eff[k, i] > 0:
                    incoming_pos_sum[i] += A_pos_eff[k, i] / out_degrees_eff[k] * x[k]
        # Compute the denominator for the transition matrix calculation, which combines the incoming positive sum and the teleportation vector z
        denom = alpha * incoming_pos_sum + (1 - alpha) * z
        denom = np.where(denom == 0, 1.0, denom)

        # Compute the transition matrix T based on the current trust scores, the adjacency matrices, and the teleportation vector z
        for i in range(n):
            for j in range(n):
                # Compute the numerator terms for the transition matrix calculation
                num_term1 = 0.0
                if A_pos_eff[j, i] > 0:
                    # Compute the contribution from positive edges, scaled by the current trust score of the source node
                    num_term1 = alpha * (A_pos_eff[j, i] / out_degrees_eff[j]) * x[j]
                # Compute the contribution from the teleportation vector z, scaled by the current trust score of the source node
                num_term2 = M * (1 - alpha) * z[i] * x[j]
                # Compute the transition probability from node j to node i by normalizing the numerator terms by the denominator
                T[i, j] = (num_term1 + num_term2) / denom[i]

        # Compute the new trust scores based on the current trust scores, the transition matrix, and the teleportation vector z
        Gx = alpha * incoming_pos_sum + (1 - alpha) * z

        # Our other new addition. For certain nodes (retracted papers) we can choose to force their trust score
        # to be 0 at every iteration. The benefit is twofold. Firstly we end up with a zero trust score for all
        # retracted papers. Secondly, retracted papers no longer confer any trust to the papers they cite.
        x_next = D_P_beta * Gx
        if forced_zero_mask.any():
            x_next[forced_zero_mask] = 0.0
        sum_x_next = np.sum(x_next)
        if sum_x_next > 0:
            x_next /= sum_x_next
        else:
            x_next = np.ones(n) / n

        # Compute the new P matrix based on the current transition matrix T and the previous P matrix
        P_tilde = T @ P

        # Compute the next P matrix based on the current P_tilde and the adjacency matrix for negative edges
        P_next = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if A_neg[i, j] > 0:
                    P_next[i, j] = 1.0
                elif i == j:
                    P_next[i, j] = 0.0
                else:
                    P_next[i, j] = P_tilde[i, j]

        # Update the diagonal of P_tilde for the next iteration
        P_tilde_diag = np.diag(P_tilde).copy()

        # Update the P matrix and the trust scores for the next iteration
        P = P_next
        x = x_next

        # Check for convergence by comparing the maximum absolute difference between the current and previous trust scores
        if np.max(np.abs(x - x_prev)) < tol:
            break

    return {nodes[i]: float(x[i]) for i in range(n)}, P
