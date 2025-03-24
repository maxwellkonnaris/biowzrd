#!/bin/bash

##################################
# 1. Configuration and Checks    #
##################################

# Repository Remote URL
# Change this repository to whatever repo you'd like
REMOTE_URL="git@github.com:Silverman-Lab/data_repository.git" 

# Ensure Git and Git LFS are installed
if ! command -v git &> /dev/null || ! command -v git-lfs &> /dev/null; then
    echo "Error: Git and/or Git LFS is not installed. Please install them and try again."
    exit 1
fi

# Ensure the current directory is a Git repository
if [ ! -d ".git" ]; then
    echo "Error: Current directory is not a Git repository. Please cd into the repository's root and try again."
    exit 1
fi

# Verify or Set the Remote Repository
CURRENT_REMOTE=$(git remote get-url origin 2>/dev/null || echo "")
if [ "$CURRENT_REMOTE" != "$REMOTE_URL" ]; then
    echo "Setting the correct remote repository..."
    git remote remove origin 2>/dev/null || echo "No existing remote to remove."
    git remote add origin "$REMOTE_URL"
fi

echo "Remote repository is set to:"
git remote -v

##############################################
# 2. Helper Function: Colored Progress Bar   #
##############################################

# show_progress "Download" $PID  (or "Upload", etc.)
# This function shows a simple progress bar while a background process ($PID) runs.
show_progress() {
    local action="$1"
    local pid="$2"
    local i=0
    local -i total=50  # length of the bar

    echo
    echo -ne "\033[93m$action in progress: ["
    # Initialize the bar with spaces
    for ((j=0; j<total; j++)); do echo -n " "; done
    echo -n "] 0%"
    echo -ne "\033[0m"

    local pos=0
    while kill -0 $pid 2>/dev/null; do
        ((pos++))
        # Calculate approximate percentage
        local percent=$((pos * 100 / total))

        # Move cursor back to start of line (carriage return) and redraw
        echo -ne "\r\033[93m$action in progress: ["
        for ((j=0; j<pos; j++)); do echo -n "#"; done
        for ((j=pos; j<total; j++)); do echo -n " "; done
        echo -n "] $percent%"
        echo -ne "\033[0m"

        sleep 0.1

        # If we reached 100% but the process is still running, we just stay at 100%
        if [ $pos -ge $total ]; then
            pos=$total
        fi
    done

    # Once the process completes, fill the bar
    echo -ne "\r\033[92m$action complete:    ["
    for ((j=0; j<total; j++)); do echo -n "#"; done
    echo -n "] 100%"
    echo -e "\033[0m\n"
}

#################################################
# 3. Fetch, Show Remote Changes, Rebase/Pull    #
#################################################

echo
echo "Stashing any local changes..."
git stash save "Backup before pull" &> /dev/null || echo "No changes to stash."

echo "Fetching latest remote data..."
# Run 'git fetch' in background to display progress bar
git fetch origin main --progress &> /dev/null &
fetch_pid=$!
show_progress "Download" $fetch_pid
wait $fetch_pid

echo "=== Remote changes (to be downloaded) ==="
# Compare HEAD..origin/main to see what changed on remote
git diff --name-status HEAD..origin/main
echo "========================================="
echo

# Check commits
LOCAL=$(git rev-parse HEAD)
REMOTE=$(git rev-parse origin/main)
BASE=$(git merge-base HEAD origin/main)

if [ "$LOCAL" = "$REMOTE" ]; then
    echo "Local branch is already up to date with remote."
elif [ "$LOCAL" = "$BASE" ]; then
    echo "Local branch is behind remote. Pulling updates..."
    git pull --rebase origin main
elif [ "$REMOTE" = "$BASE" ]; then
    echo "Local branch is ahead of remote. Ready to push."
else
    echo "Branches have diverged. Attempting to resolve with rebase..."
    if ! git pull --rebase origin main; then
        echo "Manual resolution required. Please fix conflicts."
        exit 1
    fi
fi

echo
echo "Applying any stashed changes..."
git stash pop &> /dev/null || echo "No stash to apply."

###################################
# 4. Process Repository Files     #
###################################

echo
echo "Processing repository files (zipping)..."
BASE_DIR=$(pwd)

# Find all subdirectories excluding .git
find . -type d -not -path "./.git*" | while read -r folder; do
    # Skip directories containing .gitattributes
    if [ -f "$folder/.gitattributes" ]; then
        echo "Skipping folder $folder as .gitattributes is present."
        continue
    fi

    echo "Processing folder: $folder"
    cd "$folder" || continue

    # Zip files excluding README variants, parse.R, *.zip, and *.gz
    for file in *; do
        if [ -f "$file" ]; then
            case "$file" in
                [Rr][Ee][Aa][Dd][Mm][Ee]*|parse.R|*.zip|*.gz)
                    # Skip these files
                    ;;
                *)
                    echo "Zipping $file in $folder"
                    zip -q "${file}.zip" "$file"
                    rm -f "$file"
                    ;;
            esac
        fi
    done

    cd "$BASE_DIR" || exit
done

###################################
# 5. Git LFS, Stage, Commit       #
###################################

echo
echo "Tracking .zip files with Git LFS..."
git lfs track "*.zip"

echo "Staging changes..."
git add .

# Check if there are new changes to commit
if git diff-index --quiet HEAD; then
    echo "No new changes to commit."
else
    CURRENT_DATE=$(date +"%Y-%m-%d %H:%M:%S")
    GITHUB_USERNAME=$(git config user.name)

    echo "Enter a commit message: "
    read -r COMMIT_MESSAGE

    FINAL_MESSAGE="[$CURRENT_DATE] [$GITHUB_USERNAME] $COMMIT_MESSAGE"
    git commit -m "$FINAL_MESSAGE"
fi

###################################
# 6. Show Local Changes, Push     #
###################################

# Show what will be uploaded:
# Compare origin/main..HEAD to see local changes not on remote
echo
echo "=== Local changes (to be uploaded) ==="
git diff --name-status origin/main..HEAD
echo "======================================="
echo

echo "Pushing changes to the repository..."
git push origin main --progress &> /dev/null &
push_pid=$!
show_progress "Upload" $push_pid
wait $push_pid

# If the push failed for some reason
if [ $? -ne 0 ]; then
    echo "Push failed. Attempting to rebase and push again..."
    git pull --rebase origin main && git push origin main || {
        echo "Push failed again. Manual intervention required."
        exit 1
    }
fi

echo "All done!"

