#!/bin/bash
# utils.sh -- shared utility functions for scripts

load_dependency() {
        local TOOL=$1

        if ! command -v "$TOOL" &> /dev/null; then
                echo "${TOOL} not found in PATH, attempting module load."
                if ! module load "$TOOL" 2>/dev/null; then
                        echo "ERROR: could not load ${TOOL} via module load and tool not in path."
                        echo "Please ensure ${TOOL} is installed or available as a module."
                        return 1
                fi
                if ! command -v "$TOOL" &> /dev/null; then
                        echo "ERROR: module load ${TOOL} succeeded but {$TOOL} still not found in path."
                        return 1
                fi
        else
                echo "${TOOL} found in PATH: $(command -v "$TOOL")"
        fi
        echo "${TOOL} loaded: $("$TOOL" --version 2>&1 | head -1)"

