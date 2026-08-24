#!/usr/bin/env bash
set -Eeuo pipefail

# Always run from the repository root.
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

branch="$(git branch --show-current)"

if [[ -z "$branch" ]]; then
    echo "ERROR: detached HEAD. Check out a branch before syncing."
    exit 1
fi

echo ">>> Repository: $repo_root"
echo ">>> Branch: $branch"

# Update from GitHub before creating a new commit.
echo ">>> Fetching remote changes"
git fetch origin "$branch"

echo ">>> Integrating remote changes"
git pull --rebase --autostash origin "$branch"

echo ">>> Staging tracked changes"

# Stage modifications and deletions of files already tracked by Git.
git add --update

# Explicitly stage new files only from approved locations.
allowed_paths=(
    ".gitignore"
    "README.md"
    "AGENTIC.md"
    "renv.lock"
    "renv/activate.R"
    "renv/settings.json"
    "config"
    "notes"
    "scripts"
)

for path in "${allowed_paths[@]}"; do
    if [[ -e "$path" ]]; then
        git add -- "$path"
    fi
done

# Add root-level RStudio project files, if present.
while IFS= read -r -d '' file; do
    git add -- "$file"
done < <(
    find . \
        -maxdepth 1 \
        -type f \
        -name '*.Rproj' \
        -print0
)

# Stop when there is nothing staged.
if git diff --cached --quiet; then
    echo ">>> Nothing new to commit"
    echo ">>> Synchronization completed"
    exit 0
fi

echo ">>> Changes being committed:"
git --no-pager diff --cached --stat

# Reject individual staged files larger than 10 MB.
max_size=$((10 * 1024 * 1024))

while IFS= read -r -d '' file; do
    # Deleted files no longer exist locally, so skip the size check.
    [[ -f "$file" ]] || continue

    size="$(stat -c '%s' -- "$file")"

    if (( size > max_size )); then
        size_mb="$(awk -v bytes="$size" \
            'BEGIN { printf "%.1f", bytes / 1024 / 1024 }')"

        echo
        echo "ERROR: staged file exceeds the 10 MB limit:"
        echo "  $file (${size_mb} MB)"
        echo
        echo "To unstage it:"
        printf '  git restore --staged -- %q\n' "$file"
        exit 1
    fi
done < <(
    git diff --cached \
        --name-only \
        --diff-filter=ACMR \
        -z
)

commit_message="Auto-sync $(date '+%Y-%m-%d %H:%M:%S')"

echo ">>> Creating commit"
git commit -m "$commit_message"

echo ">>> Pushing"
git push origin "$branch"

echo ">>> Synchronization completed"
