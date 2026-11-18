## Each release

- [ ] Change version number in pyinseq/VERSION
- [ ] Update CHANGELOG.md, prepare text for GitHub release
- [ ] Commit changes and ensure that the commit passes the GitHub Actions CI
- [ ] Create a GitHub release with a new tag matching pyinseq/VERSION (e.g. `0.4.0`),
      using the CHANGELOG text
- [ ] Check that the **Release** workflow in the Actions tab succeeded and the new
      version appears at https://pypi.org/project/pyinseq/

Publishing the GitHub release runs `.github/workflows/release.yml`, which builds the
package and uploads it to PyPI using Trusted Publishing; no passwords or API tokens
are needed. The workflow stops without uploading if the tag doesn't match
pyinseq/VERSION.

## One-time setup (already done once it works)

- [ ] On PyPI, as an owner of the `pyinseq` project: **Manage project → Settings →
      Publishing → Add a new publisher → GitHub**, with
      - Owner: `mjmlab`
      - Repository name: `pyinseq`
      - Workflow name: `release.yml`
      - Environment name: `pypi`
- [ ] On GitHub, in mjmlab/pyinseq **Settings → Environments**, create an environment
      named `pypi` (optionally add required reviewers so each upload needs approval)
