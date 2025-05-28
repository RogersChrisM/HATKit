import subprocess

def waitForProcess(pid, script_path="../shell/waitForProcessThenDo.sh"):
    """
    Calls shell script waitForProcessThenDo.sh with the given PID. If no PID given, assumes last created PID is PID.
    Stops calling script until process is completed.
    
    Params:
        pid (int or str): Optional process ID to wait for
        script_path (str): Path to the waitForProcessThenDo.sh script.
    """
    try:
        completed_process = subprocess.run([script_path, str(pid)], check=True, text=True, capture_output=True)
        if completed_process.stderr:
            print("ERRORS ENCOUNTERED:", completed_process.stderr)
    except subprocess.CalledProcessError as e:
        print(f"An error occured while monitoring PID {pid}: {e}")


if __name__ == '__main__':
    exit(1)
