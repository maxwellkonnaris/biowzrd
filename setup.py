import os
import stat
import argparse

def is_github_repo(repo_path):
    """Check if the given directory is a GitHub repository by looking for a .git folder."""
    return os.path.isdir(os.path.join(repo_path, ".git"))

def set_executable_permissions(repo_path):
    """Recursively finds and sets executable permissions on script files inside a GitHub repository."""
    script_extensions = {'.sh', '.py', '.pl', '.rb', '.js'}

    if not is_github_repo(repo_path):
        print(f"❌ Error: {repo_path} is not a valid GitHub repository (missing .git folder).")
        return

    print(f"🔍 Scanning GitHub repository: {repo_path}")

    for root, _, files in os.walk(repo_path):
        # Skip hidden folders (e.g., .git, .github, .venv)
        if any(folder.startswith('.') for folder in root.split(os.sep)):
            continue

        for file in files:
            if any(file.endswith(ext) for ext in script_extensions):
                file_path = os.path.join(root, file)
                os.chmod(file_path, os.stat(file_path).st_mode | stat.S_IEXEC)
                print(f"✅ Set executable: {file_path}")

    print("🎉 Permissions set for all scripts inside the GitHub repository.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Set executable permissions for all scripts in a downloaded GitHub repository.")
    parser.add_argument("repo_path", nargs="?", default=os.getcwd(), help="Path to the GitHub repository (default: current directory).")
    args = parser.parse_args()

    set_executable_permissions(args.repo_path)
