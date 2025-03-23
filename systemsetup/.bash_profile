# Source .bashrc
if [ -f ~/.bashrc ]; then
    source ~/.bashrc
fi

# Prevent running inside an already allocated interactive job
if [[ -n "$SLURM_JOB_ID" || -n "$PBS_JOBID" ]]; then
    echo "Already inside an interactive job session. Skipping job request."
    
    # Print a welcome message with useful system info
    echo -e "\n\033[1;34m-----------------------------------\033[0m"
    echo "Welcome, $(whoami)! Today is $(date)."
    echo -e "\033[1;34m-----------------------------------\033[0m"
    echo "Hostname: $(hostname)"
    echo "System Uptime: $(uptime -p)"
    echo "Current Load: $(uptime | awk -F'load average:' '{ print $2 }')"
    echo "Memory Usage: $(free -h | awk '/^Mem/ {print $3 " / " $2}')"
    echo "Disk Usage: $(df -h ~ | awk 'NR==2 {print $3 " / " $2 " used"}')"
    echo "Last Login: $(last -n 2 $USER | awk 'NR==2 {print $4, $5, $6, $7}')"
    echo -e "\033[1;34m-----------------------------------\033[0m"
       
    # Print available aliases
    echo -e "\n\033[1;34m-----------------------------------\033[0m"
    echo -e "\033[1;34mAvailable Aliases:\033[0m"
    echo -e "\033[1;34m-----------------------------------\033[0m"
    alias
    echo -e "\033[1;34m-----------------------------------\033[0m"
    
    if [[ $- == *i* ]]; then
        source /etc/profile.d/slurm.sh 2>/dev/null  # Ensure SLURM environment is loaded

        echo -e "\n\033[1;34m-----------------------------------\033[0m"
        echo -e "\033[1;34mSlurm Allocations:\033[0m"
        echo -e "\033[1;34m-----------------------------------\033[0m"
        sacctmgr show associations where user=$USER format=Account,Partition,DefaultQOS,GrpTres
        echo -e "\033[1;34m-----------------------------------\033[0m"

        echo -e "\n\033[1;34m-----------------------------------\033[0m"
        echo -e "\033[1;34mAvailable Slurm Partitions:\033[0m"
        echo -e "\033[1;34m-----------------------------------\033[0m"
        sinfo -o "%P %a %l %D %c %m %G"
        echo -e "\033[1;34m-----------------------------------\033[0m"

    fi

    # List Conda environments
    echo -e "\n\033[1;34m-----------------------------------\033[0m"
    echo -e "\033[1;34mAvailable Conda Environments:\033[0m"
    echo -e "\033[1;34m-----------------------------------\033[0m"
    conda env list 2>/dev/null || echo "Conda is not installed or not in PATH"
    echo -e "\033[1;34m-----------------------------------\033[0m"
    
    # Restore last working directory
    if [ -f ~/.last_pwd ]; then
        cd "$(cat ~/.last_pwd)" || echo "Could not restore directory"
    fi
    
    echo -e "\n\033[1;34m-----------------------------------\033[0m"
    echo -e "\033[1;34mCurrent Running Jobs (with Resources):\033[0m"
    echo -e "\033[1;34m-----------------------------------\033[0m"
    if command -v squeue &>/dev/null; then
        echo -e "JOBID              PARTITION     NAME     USER     ST    TIME      NODES NODELIST(REASON)"
        squeue -u "$USER" --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %.12R" | grep --color=auto " R " || echo -e "\033[1;33mNo running jobs.\033[0m"
    else
        echo -e "\033[1;31msqueue command not found.\033[0m"
    fi
    module load parallel
    
    return 0  # Use 'return' instead of 'exit' since this is a profile script
fi

# Ask if the user wants to request a new interactive job
read -p "Do you want to request a new interactive job? (y/n) " choice
case "$choice" in 
  y|Y ) echo "Proceeding with interactive job request...";;
  n|N ) echo "Skipping interactive job request."; return 0;;
  * ) echo "Invalid input. Exiting."; return 1;;
esac

# Check if 'requestjob' alias exists
if alias requestjob &>/dev/null; then
    # Kill the existing interactive job (if any)
    if command -v squeue &>/dev/null; then
        job_id=$(squeue -u "$USER" --noheader --format="%A %j %T" | awk '$3=="RUNNING" && $2 ~ /salloc|interactive|bash/ {print $1; exit}')
        if [[ -n "$job_id" ]]; then
            echo "Killing existing interactive job: $job_id"
            scancel "$job_id"
            sleep 2  # Allow time for cleanup
        fi
    elif command -v qstat &>/dev/null; then
        job_id=$(qstat -u "$USER" | awk '$5=="R" && $3 ~ /INTERACTIVE|bash/ {print $1; exit}')
        if [[ -n "$job_id" ]]; then
            echo "Killing existing interactive job: $job_id"
            qdel "$job_id"
            sleep 2
        fi
    fi

    # Request a new job
    echo "Requesting a new interactive job..."
    requestjob
else
    echo "Alias 'requestjob' not found!"
    return 1
fi

# Print a welcome message with useful system info
echo -e "\n\033[1;34m-----------------------------------\033[0m"
echo "Welcome, $(whoami)! Today is $(date)."
echo -e "\033[1;34m-----------------------------------\033[0m"
echo "Hostname: $(hostname)"
echo "System Uptime: $(uptime -p)"
echo "Current Load: $(uptime | awk -F'load average:' '{ print $2 }')"
echo "Memory Usage: $(free -h | awk '/^Mem/ {print $3 " / " $2}')"
echo "Disk Usage: $(df -h ~ | awk 'NR==2 {print $3 " / " $2 " used"}')"
echo "Last Login: $(last -n 2 $USER | awk 'NR==2 {print $4, $5, $6, $7}')"
echo -e "\033[1;34m-----------------------------------\033[0m"

# Print available aliases
echo -e "\n\033[1;34m-----------------------------------\033[0m"
echo -e "\033[1;34mAvailable Aliases:\033[0m"
echo -e "\033[1;34m-----------------------------------\033[0m"
alias
echo -e "\033[1;34m-----------------------------------\033[0m"

# List Conda environments
echo -e "\n\033[1;34m-----------------------------------\033[0m"
echo -e "\033[1;34mAvailable Conda Environments:\033[0m"
echo -e "\033[1;34m-----------------------------------\033[0m"
conda env list 2>/dev/null || echo "Conda is not installed or not in PATH"
echo -e "\033[1;34m-----------------------------------\033[0m"

# Restore last working directory
if [ -f ~/.last_pwd ]; then
    cd "$(cat ~/.last_pwd)" || echo "Could not restore directory"
fi

echo -e "\n\033[1;34m-----------------------------------\033[0m"
echo -e "\033[1;34mCurrent Running Jobs (with Resources):\033[0m"
echo -e "\033[1;34m-----------------------------------\033[0m"
if command -v squeue &>/dev/null; then
    echo -e "JOBID              PARTITION     NAME     USER     ST    TIME      NODES NODELIST(REASON)"
    squeue -u "$USER" --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %.12R" | grep --color=auto " R " || echo -e "\033[1;33mNo running jobs.\033[0m"
else
    echo -e "\033[1;31msqueue command not found.\033[0m"
fi
module load parallel
