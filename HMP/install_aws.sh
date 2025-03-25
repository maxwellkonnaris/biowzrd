#!/bin/bash

# Set install directories
AWS_DIR="$HOME/aws-cli"
BIN_DIR="$HOME/bin"

# Create necessary directories
mkdir -p "$AWS_DIR" "$BIN_DIR"

# Download AWS CLI
echo "Downloading AWS CLI..."
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"

# Extract AWS CLI installer
echo "Extracting AWS CLI..."
unzip -q awscliv2.zip

# Install AWS CLI to home directory
echo "Installing AWS CLI in $AWS_DIR..."
./aws/install --install-dir "$AWS_DIR" --bin-dir "$BIN_DIR"

# Remove installation files
rm -rf aws awscliv2.zip

# Add AWS CLI to PATH if not already added
if ! grep -q "$BIN_DIR" ~/.bashrc; then
    echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
    echo "Added AWS CLI to PATH in ~/.bashrc"
fi

# Apply changes to current session
export PATH="$BIN_DIR:$PATH"

# Verify installation
echo "Verifying AWS CLI installation..."
aws --version

# (Optional) Install jq (JSON processor) if needed
echo "Installing jq..."
mkdir -p "$HOME/local/bin"
curl -L -o "$HOME/local/bin/jq" "https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64"
chmod +x "$HOME/local/bin/jq"

# Add jq to PATH if not already added
if ! grep -q "$HOME/local/bin" ~/.bashrc; then
    echo 'export PATH=$HOME/local/bin:$PATH' >> ~/.bashrc
    echo "Added jq to PATH in ~/.bashrc"
fi

# Apply changes to current session
export PATH="$HOME/local/bin:$PATH"

# Verify jq installation
echo "Verifying jq installation..."
jq --version

echo "Installation complete! Please restart your terminal or run: source ~/.bashrc"
