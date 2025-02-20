#!/bin/bash
pipx-global () {
    local cmd

    cmd="$1"
    shift

    if ! python3 -m pipx --help &>/dev/null; then
        >&2 echo "Please install pipx first."
        return 1
    fi

    if pipx environment | grep -q '^PIPX_GLOBAL_HOME='; then
        python3 -m pipx "${cmd}" --global "$@" || return
    else
        PIPX_HOME="${PIPX_GLOBAL_HOME:-/opt/pipx}" \
        PIPX_BIN_DIR="${PIPX_GLOBAL_BIN_DIR:-/usr/local/bin}" \
        PIPX_MAN_DIR="${PIPX_GLOBAL_MAN_DIR:-/usr/local/share/man}" \
        python3 -m pipx "${cmd}" "$@" || return
    fi
}

pipx-global install cutadapt