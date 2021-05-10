if(BUNDLE)
  add_custom_target(html ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/Gemfile ${CMAKE_CURRENT_BINARY_DIR}/Gemfile
      COMMAND bundle install --path=vendor/bundle
      COMMAND bundle exec jekyll build --source=${CMAKE_CURRENT_SOURCE_DIR} --layouts=${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS siga
    )
endif(BUNDLE)
