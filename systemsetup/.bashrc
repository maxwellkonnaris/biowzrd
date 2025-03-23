# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific environment: Add local bin folders to PATH if not already there
if ! [[ "$PATH" =~ "$HOME/.local/bin:$HOME/bin:" ]]; then
    PATH="$HOME/.local/bin:$HOME/bin:$PATH"
fi
export PATH

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/storage/home/mak6930/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/storage/home/mak6930/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/storage/home/mak6930/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/storage/home/mak6930/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba init' !!
export MAMBA_EXE='/storage/work/mak6930/applicationstorage/bin/micromamba';
export MAMBA_ROOT_PREFIX='/storage/work/mak6930/applicationstorage/micromamba';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias micromamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

alias rm='rm -i'
alias mv='mv -i'
alias cp='cp -i'
alias ls="ls -lah --color"
alias requestjob="salloc --nodes=1 --ntasks=1 --account=one --cpus-per-task=8 --mem=150G"
alias mamba=micromamba
alias conda=micromamba

if [[ "$PBS_JOBNAME" == "jupyter-byoe-portal2" && -z "$PS1" ]]; then
    return
fi
if [[ "$PBS_JOBNAME" == "gdesktop" && -z "$PS1" ]]; then
    return
fi

export PATH=~/bin:~/edirect:$PATH
export PATH=~/bin/sratoolkit.3.0.0-ubuntu64/bin:$PATH
export PS1='\[\e]0;\w\a\]\n\[\e[01;38m\]\u \[\e[01;33m\]@\[\e[01;33m\]\h \[\e[01;37m\]\w\[\e[0m\] \[$CONDA_DEFAULT_ENV : '
export LC_ALL=C
export PATH="$HOME/bin:$PATH"
export BASH_SILENCE_DEPRECATION_WARNING=1
export PROMPT_DIRTRIM=1
export PATH=$HOME/bin:$PATH
export PATH=$HOME/local/bin:$PATH
