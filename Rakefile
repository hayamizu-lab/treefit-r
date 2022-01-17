# -*- ruby -*-

@version = nil
def version
  return @version if @version
  File.open("DESCRIPTION") do |description|
    description.each_line do |line|
      case line
      when /^Version: (.+)$/
        @version = $1
        return @version
      end
    end
  end
  nil
end

def development_version?
  version.end_with?(".9000")
end

def doc_dir
  if development_version?
    "docs/dev"
  else
    "docs/#{version}"
  end
end

def generate_document(dir)
  rm_rf(dir)
  sh("Rscript",
     "-e",
     "pkgdown::build_site(" +
     "run_dont_run=TRUE, " +
     "new_process=FALSE, " +
     "override=list(destination=\"#{dir}\"), new_process=FALSE)")
end

namespace :doc do
  desc "Generate document"
  task :generate do
    generate_document(doc_dir)
  end
end

namespace :check do
  desc "Check on win-builder.r-project.org"
  task :win_builder do
    sh("Rscript", "-e", "devtools::check_win_devel()")
    sh("Rscript", "-e", "devtools::check_win_release()")
  end
end

desc "Release"
task :release do
  sh("Rscript", "-e", "devtools::check()")
  sh("Rscript", "-e", "devtools::spell_check()")
  sh("R", "CMD", "build", ".")
  sh("R", "CMD", "check", "treefit_#{version}.tar.gz", "--as-cran")
  generate_document(doc_dir)
  sh("git", "add", doc_dir)
end

namespace :cran do
  desc "Submit to CRAN"
  task :submit do
    puts("Enter devtools::release()")
    sh("R")
  end
end

desc "Tag"
task :tag do
  sh("git", "tag",
     "-a", version,
     "-m", "Treefit #{version}!!!")
  sh("git", "push", "--tags")
end
