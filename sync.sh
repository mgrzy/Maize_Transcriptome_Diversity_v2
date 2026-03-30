# Synchronize local worktree with the remote `main` branch.

echo ">>> Pull latest changes"
# Grab upstream commits before staging local edits to minimise conflicts.
git pull origin main

echo ">>> Add and commit all changes"
# Stage every tracked/untracked change so nothing is left out of the sync.
git add -A
# Autogenerate a timestamped commit message; ignore failure when there is nothing to commit.
git commit -m "Auto-sync $(date +'%Y-%m-%d %H:%M:%S')" || true

echo ">>> Push to GitHub"
# Push local commits; emit a hint if the push fails (e.g., offline).
git push origin main || echo "Push failed, check connection"
