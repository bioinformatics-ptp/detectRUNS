The README.md is aimed at new users. The NEWS.md is aimed at existing users: it should list all the API changes in each release. There are a number of formats you can use for package news, but I recommend NEWS.md. It’s not supported by CRAN (so you’ll need to run devtools::use_build_ignore("NEWS.md")), but it’s well supported by GitHub and is easy to re-purpose for other formats.

Organise your NEWS.md as follows:

Use a top-level heading for each version: e.g. # mypackage 1.0. The most recent version should go at the top.

Each change should be included in a bulleted list. If you have a lot of changes you might want to break them up using subheadings, ## Major changes, ## Bug fixes etc. I usually stick with a simple list until just before releasing the package when I’ll reorganise into sections, if needed. It’s hard to know in advance exactly what sections you’ll need.

If an item is related to an issue in GitHub, include the issue number in parentheses, e.g. (#10). If an item is related to a pull request, include the pull request number and the author, e.g. (#101, @hadley). Doing this makes it easy to navigate to the relevant issues on GitHub.

The main challenge with NEWS.md is getting into the habit of noting a change as you make a change.

More details [here](http://r-pkgs.had.co.nz/release.html#important-files)
