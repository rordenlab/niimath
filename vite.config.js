// vite.config.js
const path = require('path')
const { defineConfig } = require('vite')

module.exports = defineConfig({
  build: {
    lib: {
      entry: path.resolve(__dirname, 'src/index.js'),
      name: 'niimath-js',
      fileName: (format) => `niimath-js.${format}.js`
    },
	}
})
