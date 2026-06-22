#!/bin/sh

fail=0 
hashes=$(git log --no-merges --format="%h" $1..$2)

for hash in $hashes; do
    # Extract the commit log
    git log -1 --format='%B' "$hash" | {
        status=0

        IFS= read -r subject
        IFS= read -r sep || sep=

        # Verify that the subject does not exceed 50 characters.
        if [ "${#subject}" -gt 50 ]; then
            echo "${hash}: Subject length ${#subject} exceeds 50 characters."
            echo "> ${subject}"
            printf '> %50s^\n' ''
            status=1
        fi

        # Verify that each line of the body does not exceed 72 characters.
        has_body=0
        line_number=0

        while IFS= read -r line
        do
            has_body=1

            line_number=$((line_number + 1))

            if [ "${#line}" -gt 72 ]; then
                echo "${hash}: Body line ${line_number} exceeds 72 characters."
                echo "> ${line}"
                printf '> %72s^\n' ''
                status=1
            fi
        done

        # Verify that the subject-body separator is blank
        if [ "${has_body}" -eq 1 ] && [ -n "${sep}" ]; then
            echo "${hash}: Missing blank line between subject and body."
            status=1
        fi

        exit "${status}"
    } || fail=1
done

exit "${fail}"
