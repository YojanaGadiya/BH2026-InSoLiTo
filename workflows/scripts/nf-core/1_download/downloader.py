
from pathlib import Path
import json
import csv
import re
import subprocess
from datetime import datetime, timezone


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parents[4]   # 1_download -> nf-core -> scripts -> workflows -> root
DATA = ROOT / "workflows" / "data" / "nf-core"

OUTPUT_ROOT = DATA / "download"              # <pipeline>/<version>/, plus manifest, snapshot, logs
CATALOGUE_PATH = OUTPUT_ROOT / "nfcore_pipelines_catalogue.json"

# None = all active pipelines in the catalogue.
# Alternatively: {"rnaseq", "sarek", "ampliseq"}
PIPELINE_NAMES = {"rnaseq", "sarek"}

# Exclude archived pipelines by default.
INCLUDE_ARCHIVED = True

# First run: True only prints the download plan.
# After checking it, change to False.
DRY_RUN = False

GIT_TIMEOUT = 900


# ============================================================
# SPARSE CHECKOUT RULES
# ============================================================

SPARSE_PATTERNS = [
    "/main.nf",

    "/workflows/",
    "/subworkflows/",

    "/nextflow.config",
    "/nextflow_schema.json",
    "/conf/",

    "/modules.json",

    # Every Nextflow file in modules/ (older pipelines name them modules/local/<process>.nf, not main.nf),
    # plus the two metadata files.
    "/*.nf",
    "/modules/**/*.nf",
    "/modules/**/meta.yml",
    "/modules/**/environment.yml",
]


# ============================================================
# HELPERS
# ============================================================

def run_git(args, cwd=None, timeout=GIT_TIMEOUT, check=True):

    result = subprocess.run(
        ["git", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout
    )

    if check and result.returncode != 0:
        raise RuntimeError(
            f"Git command failed: git {' '.join(args)}\n"
            f"{result.stderr}"
        )

    return result


def is_true(value):

    return (
        value is True
        or str(value).strip().lower() in {"true", "1", "yes"}
    )


def clean_name(value):

    return re.sub(r"[^a-zA-Z0-9._-]+", "_", str(value))


def get_latest_release(pipeline):

    releases = pipeline.get("releases") or []

    if isinstance(releases, dict):
        releases = releases.get("releases", list(releases.values()))

    if not isinstance(releases, list):
        return None

    stable = [
        r for r in releases
        if r.get("tag_sha")
        and not is_true(r.get("prerelease", False))
        and not is_true(r.get("draft", False))
    ]

    if not stable:
        return None

    # Prefer the most recently published stable release.
    dated = [r for r in stable if r.get("published_at")]

    if dated:
        return max(
            dated,
            key=lambda r: r.get("published_at", "")
        )

    return stable[0]


def load_pipeline_catalogue():

    if not CATALOGUE_PATH.exists():
        raise FileNotFoundError(
            f"Catalogue not found: {CATALOGUE_PATH}"
        )

    with open(CATALOGUE_PATH, encoding="utf-8") as f:
        catalogue = json.load(f)

    pipelines = catalogue.get("remote_workflows")

    if not isinstance(pipelines, list):
        raise ValueError(
            "Could not find a remote_workflows list in the catalogue."
        )

    return pipelines


def create_download_plan(pipelines):

    plan = []

    for pipeline in pipelines:

        full_name = (
            pipeline.get("full_name")
            or pipeline.get("name")
        )

        if not full_name:
            continue

        if "/" not in full_name:
            full_name = f"nf-core/{full_name}"

        short_name = full_name.split("/")[-1]

        if PIPELINE_NAMES is not None:
            if short_name not in PIPELINE_NAMES:
                continue

        archived = is_true(pipeline.get("archived", False))

        if archived and not INCLUDE_ARCHIVED:
            continue

        release = get_latest_release(pipeline)

        if release is None:
            continue

        version = release.get("tag_name")
        commit_sha = release.get("tag_sha")

        if not version or not commit_sha:
            continue

        plan.append({
            "pipeline": short_name,
            "full_name": full_name,
            "version": version,
            "commit_sha": commit_sha,
            "published_at": release.get("published_at"),
            "repo_url": f"https://github.com/{full_name}.git",
            "archived": archived
        })

    return sorted(plan, key=lambda x: x["pipeline"].lower())


def count_files(path):

    if not path.exists():
        return 0

    return sum(1 for p in path.rglob("*") if p.is_file())


def write_manifest(rows):

    if not rows:
        return

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    path = OUTPUT_ROOT / "download_manifest.csv"

    fields = list(rows[0].keys())

    with open(path, "w", newline="", encoding="utf-8") as f:

        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# DOWNLOAD ONE PIPELINE
# ============================================================

def download_pipeline(item):

    name = item["pipeline"]
    version = clean_name(item["version"])

    destination = OUTPUT_ROOT / name / version

    destination.mkdir(parents=True, exist_ok=True)

    git_dir = destination / ".git"

    log_messages = []

    def git(args, check=True):

        result = run_git(
            args,
            cwd=destination,
            check=check
        )

        log_messages.append(
            f"$ git {' '.join(args)}\n"
            f"{result.stdout}\n{result.stderr}"
        )

        return result

    # --------------------------------------------------------
    # Check or initialise repository
    # --------------------------------------------------------

    if git_dir.exists():

        status = git(["status", "--porcelain"])

        if status.stdout.strip():
            raise RuntimeError(
                "Existing repository has local modifications. "
                "Not overwriting it."
            )

        remote = git(
            ["remote", "get-url", "origin"],
            check=False
        )

        if remote.returncode == 0:

            existing_url = remote.stdout.strip()

            if existing_url.rstrip(".git") != item["repo_url"].rstrip(".git"):
                raise RuntimeError(
                    f"Unexpected existing remote: {existing_url}"
                )

        else:
            git(["remote", "add", "origin", item["repo_url"]])

    else:

        if any(destination.iterdir()):
            raise RuntimeError(
                "Destination exists and is not an empty Git repository."
            )

        git(["init", "--quiet"])
        git(["remote", "add", "origin", item["repo_url"]])

    # --------------------------------------------------------
    # Configure sparse checkout BEFORE checking out files
    # --------------------------------------------------------

    git(["config", "core.sparseCheckout", "true"])
    git(["config", "core.sparseCheckoutCone", "false"])

    sparse_file = git_dir / "info" / "sparse-checkout"
    sparse_file.parent.mkdir(parents=True, exist_ok=True)

    sparse_file.write_text(
        "\n".join(SPARSE_PATTERNS) + "\n",
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # Fetch exactly one commit, without downloading full history
    # --------------------------------------------------------

    fetch = git([
        "fetch",
        "--depth=1",
        "--filter=blob:none",
        "origin",
        item["commit_sha"]
    ], check=False)

    if fetch.returncode != 0:

        # Fallback to the release tag if fetching by SHA fails.
        fetch = git([
            "fetch",
            "--depth=1",
            "--filter=blob:none",
            "origin",
            f"refs/tags/{item['version']}"
        ], check=False)

        if fetch.returncode != 0:
            raise RuntimeError(
                "Could not fetch the commit or release tag."
            )

    # --------------------------------------------------------
    # Checkout the exact commit
    # --------------------------------------------------------

    git(["checkout", "--force", "--detach", "FETCH_HEAD"])

    actual_sha = git(["rev-parse", "HEAD"]).stdout.strip()

    if actual_sha.lower() != item["commit_sha"].lower():
        raise RuntimeError(
            f"Commit mismatch: expected {item['commit_sha']}, "
            f"got {actual_sha}"
        )

    shallow = git(
        ["rev-parse", "--is-shallow-repository"]
    ).stdout.strip()

    # --------------------------------------------------------
    # Count downloaded contents
    # --------------------------------------------------------

    module_main = list(destination.glob("modules/**/main.nf"))
    module_meta = list(destination.glob("modules/**/meta.yml"))
    module_env = list(destination.glob("modules/**/environment.yml"))

    summary = {
        **item,
        "actual_commit": actual_sha,
        "shallow_repository": shallow,
        "destination": str(destination),
        "workflow_files": count_files(destination / "workflows"),
        "subworkflow_files": count_files(destination / "subworkflows"),
        "module_main_files": len(module_main),
        "module_meta_files": len(module_meta),
        "module_environment_files": len(module_env),
        "total_selected_files": (
            count_files(destination)
        ),
        "status": "downloaded"
    }

    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_path = log_dir / f"{name}_{version}.log"

    log_path.write_text(
        "\n\n".join(log_messages),
        encoding="utf-8"
    )

    return summary


# ============================================================
# MAIN
# ============================================================

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

print("Loading nf-core pipeline catalogue...")

catalogue_pipelines = load_pipeline_catalogue()

download_plan = create_download_plan(catalogue_pipelines)

print(f"\nCatalogue entries: {len(catalogue_pipelines)}")
print(f"Pipelines selected: {len(download_plan)}")

print("\nDownload plan:")
for item in download_plan:
    print(
        f"{item['full_name']:30} "
        f"{item['version']:15} "
        f"{item['commit_sha']}"
    )

# Preserve the catalogue used for this download.
with open(
    OUTPUT_ROOT / "catalogue_snapshot.json",
    "w",
    encoding="utf-8"
) as f:
    json.dump(catalogue_pipelines, f, indent=2)

if DRY_RUN:

    print("\nDRY RUN ONLY.")
    print("Check the selected pipelines above.")
    print("When ready, set DRY_RUN = False and rerun this cell.")

else:

    results = []

    for i, item in enumerate(download_plan, start=1):

        print(
            f"\n[{i}/{len(download_plan)}] "
            f"Downloading {item['full_name']} "
            f"{item['version']}"
        )

        try:
            result = download_pipeline(item)
            results.append(result)

            print(
                f"Completed | "
                f"Modules: {result['module_main_files']} | "
                f"Files: {result['total_selected_files']}"
            )

        except Exception as error:

            failed = {
                **item,
                "destination": str(
                    OUTPUT_ROOT
                    / item["pipeline"]
                    / clean_name(item["version"])
                ),
                "status": "failed",
                "error": str(error)
            }

            results.append(failed)

            print("FAILED:", error)

        # Save progress after every pipeline.
        write_manifest(results)

    # Final report
    successful = [
        r for r in results
        if r["status"] == "downloaded"
    ]

    failed = [
        r for r in results
        if r["status"] == "failed"
    ]

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "catalogue_path": str(CATALOGUE_PATH),
        "output_path": str(OUTPUT_ROOT),
        "selected_pipelines": len(download_plan),
        "downloaded": len(successful),
        "failed": len(failed),
        "total_module_main_files": sum(
            r.get("module_main_files", 0)
            for r in successful
        ),
        "total_selected_files": sum(
            r.get("total_selected_files", 0)
            for r in successful
        )
    }

    with open(
        OUTPUT_ROOT / "download_summary.json",
        "w",
        encoding="utf-8"
    ) as f:
        json.dump(report, f, indent=2)

    print("\n" + "=" * 55)
    print("DOWNLOAD COMPLETED")
    print("=" * 55)

    for key, value in report.items():
        print(f"{key}: {value}")
