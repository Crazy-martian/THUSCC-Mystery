#!/usr/bin/env bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 [--dry] <problem> <file>"
    exit 1
fi

DRY_RUN="false"
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry)
            DRY_RUN=true
            shift
            ;;
        *)
            if [ -z "$PROBLEM" ]; then
                PROBLEM="$1"
            elif [ -z "$FILE" ]; then
                FILE="$1"
            else
                echo "Error: Too many arguments."
                echo "Usage: $0 [--dry] <problem> <file>"
                exit 1
            fi
            shift
            ;;
    esac
done
if [ -z "$PROBLEM" ]; then
    echo "Error: Problem parameter is required."
    echo "Usage: $0 [--dry] <problem> <file>"
    exit 1
fi
if [ -z "$FILE" ]; then
    echo "Error: File parameter is required."
    echo "Usage: $0 [--dry] <problem> <file>"
    exit 1
fi

# Credentials: support env overrides, project or home credentials file.
# File format (recommended):
# username=your_username
# password=your_password
# If file doesn't exist or is missing fields, prompt interactively and save.

# Environment overrides (take precedence)
USERNAME="${THUSCC_USERNAME:-}"
PASSWORD="${THUSCC_PASSWORD:-}"

# Credential file selection: explicit env, project .thuscc_credentials, then $HOME
if [ -n "${THUSCC_CREDENTIALS_FILE:-}" ]; then
    CREDENTIALS_FILE="$THUSCC_CREDENTIALS_FILE"
elif [ -f "./.thuscc_credentials" ]; then
    CREDENTIALS_FILE="./.thuscc_credentials"
else
    CREDENTIALS_FILE="$HOME/.thuscc_credentials"
fi

# Read from credentials file if exists and values not provided by env
if [ -f "$CREDENTIALS_FILE" ]; then
    # Try key=value parsing
    if [ -z "$USERNAME" ]; then
        USERNAME=$(grep -E '^username[[:space:]]*=' "$CREDENTIALS_FILE" 2>/dev/null | head -n1 | cut -d'=' -f2- | sed 's/^\s*//;s/\s*$//')
    fi
    if [ -z "$PASSWORD" ]; then
        PASSWORD=$(grep -E '^password[[:space:]]*=' "$CREDENTIALS_FILE" 2>/dev/null | head -n1 | cut -d'=' -f2-)
    fi
    # Fallback: if file has no key= lines, treat first line as username and second as password
    if [ -z "$USERNAME" ] || [ -z "$PASSWORD" ]; then
        LINE1=$(sed -n '1p' "$CREDENTIALS_FILE" 2>/dev/null || true)
        LINE2=$(sed -n '2p' "$CREDENTIALS_FILE" 2>/dev/null || true)
        if [ -z "$USERNAME" ] && [ -n "$LINE1" ]; then
            USERNAME="$LINE1"
        fi
        if [ -z "$PASSWORD" ] && [ -n "$LINE2" ]; then
            PASSWORD="$LINE2"
        fi
    fi
fi

# If still missing, prompt interactively
CREATED_CRED_FILE=false
if [ -z "$USERNAME" ]; then
    read -p "Enter your submit username: " USERNAME
    CREATED_CRED_FILE=true
fi
if [ -z "$PASSWORD" ]; then
    read -sp "Enter your submit password: " PASSWORD
    echo
    CREATED_CRED_FILE=true
fi

# Save credentials if we collected them interactively and file does not exist or is not writable
if [ "$CREATED_CRED_FILE" = true ]; then
    # Ensure parent dir exists for given path
    mkdir -p "$(dirname "$CREDENTIALS_FILE")" 2>/dev/null || true
    printf "username=%s\npassword=%s\n" "$USERNAME" "$PASSWORD" > "$CREDENTIALS_FILE"
    chmod 600 "$CREDENTIALS_FILE" 2>/dev/null || true
    printf "\033[0;32mCredentials saved to file '$CREDENTIALS_FILE'\033[0m\n"
fi

if [ ! -f "$FILE" ]; then
    echo "Error: File '$FILE' does not exist."
    exit 1
fi
if ! file "$FILE" | grep -q "text"; then
    echo "Error: File '$FILE' is not a text file."
    exit 1
fi

TZ=Asia/Shanghai
SUBMIT_TIME=$(env TZ=$TZ date +"%Y-%m-%d %H:%M:%S")
echo "Submission Time ($TZ): $SUBMIT_TIME"

SERVER_URL="https://oj.cs.tsinghua.edu.cn/thuscc/api/submit"

# Build JSON payload safely (simple quoting)
PAYLOAD=$(cat <<EOF
{
  "username": "${USERNAME}",
  "password": "${PASSWORD}",
  "problem": "${PROBLEM}",
  "submitTime": "${SUBMIT_TIME}",
  "dryRun": ${DRY_RUN}
}
EOF
)

RESPONSE=$(curl -s -X POST \
    -F "file=@$FILE" \
    -F "data=$PAYLOAD" "$SERVER_URL")

MESSAGE=$(echo "$RESPONSE" | grep -oP '"message"\s*:\s*"\K[^"]+' || true)

if [ -n "$MESSAGE" ]; then
    printf "Server Message: "
    if [ "$MESSAGE" = "Validation Passed (dry run)" ] || [ "$MESSAGE" = "Submission Accepted" ]; then
        printf "\033[1;32m%s\033[0m\n" "$MESSAGE"
    else
        printf "\033[1;31m%s\033[0m\n" "$MESSAGE"
    fi

    case "$MESSAGE" in
        "Wrong Password"|"Unknown Team")
            if [ -f "$CREDENTIALS_FILE" ]; then
                rm -f "$CREDENTIALS_FILE"
                printf "\033[0;31mRemoved invalid credentials file: '%s'\033[0m\n" "$CREDENTIALS_FILE"
                printf "\033[0;31mPlease re-run the submission script to enter correct credentials\033[0m\n"
            fi
            ;;
    esac
else
    echo "Error: Failed to parse server response."
    echo "Full response: $RESPONSE"
fi