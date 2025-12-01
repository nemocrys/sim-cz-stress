import os
import subprocess
import sys
import time

def run_step(directory, script_name):
    original_cwd = os.getcwd()
    target_dir = os.path.abspath(directory)
    
    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' not found.")
        sys.exit(1)

    try:
        print(f"Changing directory to: {target_dir}")
        os.chdir(target_dir)

        print(f"Running {script_name}...")
        result = subprocess.run([sys.executable, script_name], capture_output=False)
        
        if result.returncode != 0:
            print(f"Error: {script_name} failed with return code {result.returncode}.")
            sys.exit(result.returncode)
        else:
            print(f"Successfully ran {script_name}.")

    except Exception as e:
        print(f"An unexpected error occurred while running {script_name}: {e}")
        sys.exit(1)
    finally:
        os.chdir(original_cwd)
        print(f"Restored working directory to: {original_cwd}")


def main():

    real_start = time.time()
    cpu_start = os.times() 


    # Step 1: Run 2D setup
    print("\n==============================")
    print("RUNNING 2D SIMULATION")
    print("==============================")
    run_step(os.path.join("2D", "Csi_reference_case"), "setup.py")
    # run_step(os.path.join("2D", "Csi_optimum"), "setup.py")
    

    # Step 2: Run 3D simulation
    print("\n==============================")
    print("RUNNING 3D SIMULATION")
    print("==============================")
    run_step("3D", "run.py")



    # End timers
    real_end = time.time()
    cpu_end = os.times()
    total_real = real_end - real_start
    total_cpu = (cpu_end.user + cpu_end.system) - (cpu_start.user + cpu_start.system)

    print("\n==============================")
    print("TOTAL EXECUTION TIME SUMMARY")
    print("==============================")
    print(f"Total real time: {total_real:.2f} seconds ({total_real/60:.2f} minutes)")
    print(f"Total CPU time:        {total_cpu:.2f} seconds ({total_cpu/60:.2f} minutes)")
    print("==============================\n")

if __name__ == "__main__":
    main()
