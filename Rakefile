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

def doc_dir
  "docs/#{version}"
end

desc "Release"
task :release do
  sh("Rscript", "-e", "devtools::check()")
  sh("Rscript", "-e", "devtools::spell_check()")
  sh("R", "CMD", "build", ".")
  sh("R", "CMD", "check", "treefit_#{version}.tar.gz", "--as-cran")
  rm_rf(doc_dir)
  sh("Rscript",
     "-e",
     "pkgdown::build_site(override=list(destination=\"#{doc_dir}\"))")
  sh("git", "add", doc_dir)
end

desc "Tag"
task :tag do
  sh("git", "tag", "-a", version, "-m", "Publish #{version}")
  sh("git", "push", "--tags")
end
