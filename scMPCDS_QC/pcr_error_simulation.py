import numpy as np
import matplotlib.pyplot as plt

# random seed
np.random.seed(42)

def simulate_pcr_error(cycles=30, fidelity=1e-6, initial_copies=1, target_length=1000):
    """
    Using expectation method to mimic the growth of nucleotide error rate during PCR duplication
    
    Parameter:
    cycles: number of PCR cycles
    fidelity: high fidelity enzyme ( error_rate = 1 - fidelity)
    initial_copies: the initial template number
    target_length: the base pair lengths of target duplication fragments
    
    返回:
    results: 包含每个循环信息的列表
    """
    error_rate = 1 - fidelity
    results = []
    
    # initializing status
    total_copies = initial_copies
    total_errors = 0
    total_bases = initial_copies * target_length
    
    for cycle in range(1, cycles + 1):
        # assume that the duplication efficiency is 100% throughout the entire PCR assay.
        new_copies = total_copies
        total_copies *= 2
        
        # calculate the error number induced by nascent synthesized chains
        bases_synthesized = new_copies * target_length
        
        # use the expected error to avoid excessively large value
        expected_errors = bases_synthesized * error_rate
        
        # introduce stochastic fluctuation
        if expected_errors > 0 and expected_errors < 1e9:
            errors_this_cycle = np.random.poisson(expected_errors)
        else:
            # use the expected number if the expected error is too large
            errors_this_cycle = int(np.round(expected_errors))
        
        total_errors += errors_this_cycle
        total_bases = total_copies * target_length
        
        # calculate error ratio
        error_ratio = total_errors / total_bases if total_bases > 0 else 0
        
        results.append({
            'cycle': cycle,
            'copies': total_copies,
            'errors': total_errors,
            'total_bases': total_bases,
            'error_ratio': error_ratio
        })
    
    return results

def plot_pcr_error(results, output_file='pcr_error_plot.png'):
    """
    Visualization of error accumunition of PCR assay
    
    Parameter:
    results: Simulation results
    output_file: output results
    """
    cycles = [r['cycle'] for r in results]
    error_ratios = [r['error_ratio'] for r in results]
    errors = [r['errors'] for r in results]
    copies = [r['copies'] for r in results]
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    
    # Variation of Error ratio alongside cell cycle duplication
    ax1.plot(cycles, error_ratios, marker='o', linestyle='-', color='red')
    ax1.set_xlabel('PCR Cycle')
    ax1.set_ylabel('Error Ratio')
    ax1.set_title('Error Ratio vs PCR Cycles')
    ax1.grid(True)
    
    # Variation of Error numbers alongside cell cycle duplication
    ax2.plot(cycles, errors, marker='o', linestyle='-', color='blue')
    ax2.set_xlabel('PCR Cycle')
    ax2.set_ylabel('Total Errors')
    ax2.set_title('Total Errors vs PCR Cycles')
    ax2.grid(True)
    
    # Variation of copy numbers alongside cell cycle duplication
    ax3.plot(cycles, copies, marker='o', linestyle='-', color='green')
    ax3.set_xlabel('PCR Cycle')
    ax3.set_ylabel('Number of Copies')
    ax3.set_title('Amplicon Copies vs PCR Cycles')
    ax3.grid(True)
    ax3.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def main():
    # Error rate of High-fidelity DNA Polymerase
    high_fidelity_1 = 1 - 1e-6
    
    results1 = simulate_pcr_error(cycles=70, fidelity=high_fidelity_1)
    
    # Summary report
    final1 = results1[-1]
    print(f"After 70 PCR duplication cycles:")
    print(f"  Total copy number: {final1['copies']:,}")
    print(f"  Total nucleotide number: {final1['total_bases']:,}")
    print(f"  Total error: {final1['errors']}")
    print(f"  Error ratio: {final1['error_ratio']:.6f}")
    print(f"  Error rate: {(final1['errors'] / final1['total_bases'] * 100):.8f}%")
    
    # Viualize simulation results
    plot_pcr_error(results1, 'pcr_error_plot_1e-6_v70.pdf')
    print("Figure saved as: pcr_error_plot_1e-6_v70.pdf")

if __name__ == "__main__":
    main()
